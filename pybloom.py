#coding:utf8
"""This module implements a bloom filter probabilistic data structure and
an a Scalable Bloom Filter that grows in size as your add more items to it
without increasing the false positive error_rate.

Requires the bitarray library: http://pypi.python.org/pypi/bitarray/

    >>> from pybloom import BloomFilter
    >>> f = BloomFilter(capacity=10000, error_rate=0.001)
    >>> for i in xrange(0, f.capacity):
    ...     _ = f.add(i)
    ...
    >>> 0 in f
    True
    >>> f.capacity in f
    False
    >>> len(f) <= f.capacity
    True
    >>> abs((len(f) / float(f.capacity)) - 1.0) <= f.error_rate
    True

    >>> from pybloom import ScalableBloomFilter
    >>> sbf = ScalableBloomFilter(mode=ScalableBloomFilter.SMALL_SET_GROWTH)
    >>> count = 10000
    >>> for i in xrange(0, count):
    ...     _ = sbf.add(i)
    ...
    >>> sbf.capacity > count
    True
    >>> len(sbf) <= count
    True
    >>> abs((len(sbf) / float(count)) - 1.0) <= sbf.error_rate
    True

"""
import math
import hashlib
from struct import unpack, pack

try:
    import bitarray
except ImportError:
    raise ImportError('pybloom requires bitarray >= 0.3.4')

__version__ = '1.0.2'
__author__  = "Jay Baird <jay@mochimedia.com>, Bob Ippolito <bob@redivi.com>"

def make_hashfuncs(num_slices, num_bits):
    if num_bits >= (1 << 31):
        fmt_code, chunk_size = 'Q', 8
    elif num_bits >= (1 << 15):
        fmt_code, chunk_size = 'I', 4
    else:
        fmt_code, chunk_size = 'H', 2
    total_hash_bits = 8 * num_slices * chunk_size
    if total_hash_bits > 384:
        hashfn = hashlib.sha512
    elif total_hash_bits > 256:
        hashfn = hashlib.sha384
    elif total_hash_bits > 160:
        hashfn = hashlib.sha256
    elif total_hash_bits > 128:
        hashfn = hashlib.sha1
    else:
        hashfn = hashlib.md5
    fmt = fmt_code * (hashfn().digest_size // chunk_size)
    num_salts, extra = divmod(num_slices, len(fmt))
    if extra:
        num_salts += 1
    salts = [hashfn(hashfn(pack('I', i)).digest()) for i in xrange(num_salts)]
    def _make_hashfuncs(key):
        if isinstance(key, basestring) and not isinstance(key, unicode):
            key = key.encode('utf-8')
        else:
            key = str(key)
        rval = []
        for salt in salts:
            h = salt.copy()
            h.update(key)
            rval.extend(uint % num_bits for uint in unpack(fmt, h.digest()))
        del rval[num_slices:]
        return rval
    return _make_hashfuncs


class BloomFilter(object):
    def __init__(self, capacity, error_rate=0.01):
        print "Initialising pybloom BF(capacity="+str(capacity)+", error_rate="+str(error_rate)+")"
        """Implements a space-efficient probabilistic data structure

        capacity
            this BloomFilter must be able to store at least *capacity* elements
            while maintaining no more than *error_rate* chance of false
            positives
        error_rate
            the error_rate of the filter returning false positives. This
            determines the filters capacity. Inserting more than capacity
            elements greatly increases the chance of false positives.

        >>> b = BloomFilter(capacity=100000, error_rate=0.001)
        >>> b.add("test")
        False
        >>> "test" in b
        True

        """
        if not (0 < error_rate < 1):
            raise ValueError("Error_Rate must be between 0 and 1.")
        if not capacity > 0:
            raise ValueError("Capacity must be > 0")
        # given M = num_bits, k = num_slices, p = error_rate, n = capacity
        # solving for m = bits_per_slice
        # n ~= M * ((ln(2) ** 2) / abs(ln(P)))
        # n ~= (k * m) * ((ln(2) ** 2) / abs(ln(P)))
        # m ~= n * abs(ln(P)) / (k * (ln(2) ** 2))
        num_slices = int(math.ceil(math.log(1 / error_rate, 2)))
        # the error_rate constraint assumes a fill rate of 1/2
        # so we double the capacity to simplify the API
        bits_per_slice = int(math.ceil(
            (2 * capacity * abs(math.log(error_rate))) /
            (num_slices * (math.log(2) ** 2))))
        self.error_rate = error_rate
        self.num_slices = num_slices
        self.bits_per_slice = bits_per_slice
        self.capacity = capacity
        self.num_bits = num_slices * bits_per_slice
        self.count = 0
        #print '\n'.join('%s = %s' % tpl for tpl in sorted(self.__dict__.items()))
        self.bitarray = bitarray.bitarray(self.num_bits)
        self.bitarray.setall(False)
        self.make_hashes = make_hashfuncs(self.num_slices, self.bits_per_slice)

    def __contains__(self, key):
        """Tests a key's membership in this bloom filter.

        >>> b = BloomFilter(capacity=100)
        >>> b.add("hello")
        False
        >>> "hello" in b
        True

        """
        bits_per_slice = self.bits_per_slice
        bitarray = self.bitarray
        if not isinstance(key, list):
            hashes = self.make_hashes(key)
        else:
            hashes = key
        offset = 0
        for k in hashes:
            if not bitarray[offset + k]:
                return False
            offset += bits_per_slice
        return True

    def __len__(self):
        """Return the number of keys stored by this bloom filter."""
        return self.count

    def add(self, key, skip_check=False):
        """ Adds a key to this bloom filter. If the key already exists in this
        filter it will return True. Otherwise False.

        >>> b = BloomFilter(capacity=100)
        >>> b.add("hello")
        False
        >>> b.add("hello")
        True

        """
        bitarray = self.bitarray
        bits_per_slice = self.bits_per_slice
        hashes = self.make_hashes(key)
        if not skip_check and hashes in self:
            return True
        if self.count > self.capacity:
            print self
            self.print_bitarray()
            raise IndexError("BloomFilter is at capacity")
        offset = 0
        for k in hashes:
            self.bitarray[offset + k] = True
            offset += bits_per_slice
        self.count += 1
        return False


    #-------------------------Föll frá mér. Þetta eru einu breytingarnar á klasanum------------------------
    def computeRatio(self):
        cZero = 0
        cOne = self.bitarray.count()
        cZero = self.bitarray.length() - cOne
        assert (cOne+cZero)==self.bitarray.length()

        ratio = float(cOne)/(cZero+cOne)
        if ratio >= 0.5:
            assert cOne>=cZero
        else:
            assert cOne<cZero

        return cZero,cOne,ratio

    def print_bitarray(self):
        cZero,cOne,ratio = self.computeRatio()
        print "\nPrinting pybloom.bitArray:" \
        +"\n  Number of ones:  "+str(cOne) \
        +"\n  Number of zeros: "+str(cZero) \
        +"\n  Ratio:           "+str(ratio)

    def __str__(self):
        print "\nPrinting the pybloom BF (except ratio):"
        return \
           "  p:          "+str(self.error_rate) \
        +"\n  m:          "+str(self.bitarray.length()) \
        +"\nAdded "+str(self.count)+" kmers\n"





if __name__ == "__main__":
    import doctest
    doctest.testmod()
