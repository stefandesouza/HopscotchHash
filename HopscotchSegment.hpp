#pragma once

#include <iostream>
#include <utility>
#include <limits>
#include <cstring>
#include <climits>
#include <memory>
#include <thread>
#include <mutex>
#include <cassert>
#include <atomic>
#include <thread>
#include <stdexcept>

#define HYPRE_XXH_rotl32(x,r) ((x << r) | (x >> (32 - r)))
#define HYPRE_XXH_PRIME32_1   2654435761U
#define HYPRE_XXH_PRIME32_2   2246822519U
#define HYPRE_XXH_PRIME32_3   3266489917U
#define HYPRE_XXH_PRIME32_4    668265263U
#define HYPRE_XXH_PRIME32_5    374761393U

namespace GaspiLS {


template <typename IT, typename S>

class HopscotchSegment{

private:

  static const short _NULL_DELTA = SHRT_MIN;
  static const IT    _NULL_KEY   = -1;
  static const S     _NULL_VALUE = 0;

  struct Bucket
  {
    short _first_delta;     // first location that hashes to the home bucket
    short _next_delta;      // location of next element in the chain
    std::atomic<IT> _key;   // entry to be inserted
    S _value;               // value associated with the key
    IT _hash; 		    // bucket hash

    void init(const IT &hashValue)
    {
      _first_delta = _NULL_DELTA;
      _next_delta = _NULL_DELTA;
      _key = _NULL_KEY;
      _value = _NULL_VALUE;
      _hash = hashValue;
    }
  };

  struct Segment
  {
    IT _timestamp;
    std::mutex _mtx; // mutex to lock the segment

    void init()
    {
	_timestamp = 0;
    }
  };

  std::unique_ptr<Bucket[]>  _table;
  std::unique_ptr<Segment[]> _segments;

  IT _size;
  IT _hop_range;
  IT _concurrency_level;

  void
  initializeTable();

  inline IT
  hashFunction(const IT &key) const
  {
    
    int h32 = HYPRE_XXH_PRIME32_5 + sizeof(key);

    // 1665863975 is added to input so that
    // only -1073741824 gives HYPRE_HOPSCOTCH_HASH_EMPTY.
    // Hence, we're fine as long as key is non-negative.
    h32 += (key + 1665863975)*HYPRE_XXH_PRIME32_3;
    h32 = HYPRE_XXH_rotl32(h32, 17)*HYPRE_XXH_PRIME32_4;

    h32 ^= h32 >> 15;
    h32 *= HYPRE_XXH_PRIME32_2;
    h32 ^= h32 >> 13;
    h32 *= HYPRE_XXH_PRIME32_3;
    h32 ^= h32 >> 16;

    //hypre_assert(HYPRE_HOPSCOTCH_HASH_EMPTY != h32);

    //return h32;

    return (h32 & (_size - 1));
    
    //return (key * 2) % _size;

  }

  inline IT
  segmentFunction(const IT &hash) const
  {
    return (hash % (_concurrency_level));
    //return hash % _concurrency_level;
  }

  inline IT
  powerOfTwo(const IT &n) const
  {
    IT p = 1;
    if (n && !(n & (n - 1)))
        return n;
 
    while (p < n)
        p <<= 1;
     
    return p;
  }

  bool
  linearProbe
    ( const IT key
    , Bucket **freeBucket
    , Bucket **prevBucket
    , IT & freeDist
    );

  void
  remove_key
    ( Bucket *const from_bucket
    , Bucket *const key_bucket
    , Bucket *const prev_key_bucket
    , const IT hash );

  bool
  findCloserFreeBucket
    ( const IT &key
    , Segment * segment
    , Bucket **freeBucket
    , IT & freeDist
    );

  void
  swapPositions
    ( Segment * segment
    , Bucket **freeSpace
    , Bucket *swapSpace );

  void
  addToNeighborhood
    ( const IT key
    , Bucket *freeBucket
    , Bucket *prevBucket
    , const S value);

  void
  advance
    ( Bucket **current_bucket
    , const IT next_delta) const;

  IT
  distance
    ( const IT & hashBeg
    , const IT & hashEnd ) const;

public:
  HopscotchSegment();

  HopscotchSegment(const IT &size);

  HopscotchSegment(const IT &size, const IT &hop_range);

  HopscotchSegment(const IT &size, const IT &hop_range, const IT &concurrency_level);

  HopscotchSegment(const HopscotchSegment & other);

  ~HopscotchSegment();

  inline void
  reserve(const IT &size);

  inline void
  reserve(const IT &size, const IT &hop_range);

  inline void
  reserve(const IT &size, const IT &hop_range, const IT &concurrency_level);

  inline bool
  contains(const IT &key) const;

  inline IT
  find(const IT &key) const;

  std::pair<IT,bool>
  insert(const IT &key, const S &value);

  std::pair<IT,bool>
  set(const IT &key, const S &value);

  bool
  remove(const IT &key);

  void
  displayTable();

  void
  displayLists();
        

};

////////////////////////////////////////////////////////////////////////////////

// empty constructor
template <typename IT, typename S>
HopscotchSegment<IT, S>
    ::HopscotchSegment()
: HopscotchSegment(10)
{ }

// constructor with just the size
template <typename IT, typename S>
HopscotchSegment<IT, S>
    ::HopscotchSegment(const IT &size)
: HopscotchSegment(size,size/5)
{ }

// constructor with size and hop range
template <typename IT, typename S>
HopscotchSegment<IT, S>
    ::HopscotchSegment(const IT &size, const IT &hop_range)
: HopscotchSegment(size,hop_range,std::thread::hardware_concurrency())
{ }

// constructor with size, hop range and concurrency level
template <typename IT, typename S>
HopscotchSegment<IT, S>
    ::HopscotchSegment
      (const IT &size, const IT &hop_range, const IT &concurrency_level)
{
    _hop_range = hop_range;
    _size = powerOfTwo(size);
    _concurrency_level = powerOfTwo(concurrency_level);
    initializeTable();
}

// Destructor
template <typename IT, typename S>
HopscotchSegment<IT, S>
    ::~HopscotchSegment()
{
   _table.reset();
   _segments.reset();
}

template <typename IT, typename S>
void
HopscotchSegment<IT, S>
  ::reserve(const IT& size)
{
    reserve(size, size, 12);
}

// constructor with size, hop range and concurrency level
template <typename IT, typename S>
void
HopscotchSegment<IT, S>
    ::reserve
      (const IT &size, const IT &hop_range)
{
    reserve(size, hop_range, 12);
}

// constructor with size, hop range and concurrency level
template <typename IT, typename S>
void
HopscotchSegment<IT, S>
    ::reserve
      (const IT &size, const IT &hop_range, const IT &concurrency_level)
{
    _hop_range = hop_range;
    _size = powerOfTwo(size);
    _concurrency_level = powerOfTwo(concurrency_level);
    initializeTable();
}


// initialize table of hashBuckets
template <typename IT, typename S>
void 
HopscotchSegment<IT, S>
    ::initializeTable()
{
    _table.reset(new Bucket[_size]);
    _segments.reset(new Segment[_concurrency_level]);


    for (IT i = 0; i < this->_size; ++i)
    {
      _table[i].init(i);
    }

    for (IT i = 0; i < _concurrency_level; ++i)
    {
      _segments[i].init();
    }

}

template <typename IT, typename S>
HopscotchSegment<IT, S>
  ::HopscotchSegment(const HopscotchSegment & other)
: HopscotchSegment()
{
  _size = other._size;
  _hop_range = other._hop_range;
  _concurrency_level = other._concurrency_level;

  _table.reset(new Bucket[_size]);
  _segments.reset(new Segment[_concurrency_level]);

  for(IT i(0);i<_size;++i) {
    _table[i]._first_delta = other._table[i]._first_delta;
    _table[i]._next_delta = other._table[i]._next_delta;
    _table[i]._key.store(other._table[i]._key, std::memory_order_relaxed);
    _table[i]._value = other._table[i]._value;
    _table[i]._hash = i;
  }

}

// find an open spot in the neighborhood
template <typename IT, typename S>
bool
HopscotchSegment<IT, S>
    ::linearProbe
      ( const IT key
      , Bucket **freeBucket
      , Bucket **prevBucket
      , IT & freeDist
      )
{
  IT hash(hashFunction(key));

  *prevBucket = NULL;
  *freeBucket = &_table[hash];
  freeDist = 0;
  IT nullkey = _NULL_KEY;

  bool loop = false;

  while (!((*freeBucket)->_key == _NULL_KEY && // atomic check for free bucket's key
           (*freeBucket)->_key.compare_exchange_weak
	       ( nullkey, key, std::memory_order_release
	       , std::memory_order_relaxed))
        )
  {
    // need to keep track if a previous bucket hashes to the same value so
    // it's next delta can be updated. I am the only thread in this segment
    // i.e. no atomic check required here
    if (hashFunction((*freeBucket)->_key) == hash)
    {
      *prevBucket = *freeBucket;
    }

    ++(*freeBucket);++freeDist;

    // if size is exceeded, loop to start
    if ((*freeBucket)-_table.get() >= _size)
    {
      *freeBucket = &(_table[0]);
      loop = true;
      continue;
    }
    // no open spots available, a resize is required
    else if (loop == true && (*freeBucket)->_hash == hash)
    {
      std::cout << _size <<"\n";
      std::cout << freeDist<<"\n";
      return false;
      
    }
    
  }
  // free bucket found
  return true;
}

// remove function to set the deltas right
template <typename IT, typename S>
void
HopscotchSegment<IT, S>
  ::remove_key
    ( Bucket *const home_bucket
    , Bucket *const key_bucket
    , Bucket *const prev_key_bucket
    , const IT hash
    )
{
  key_bucket->_key.store(_NULL_KEY, std::memory_order_relaxed);
  key_bucket->_value = _NULL_VALUE;

  if (NULL == prev_key_bucket) {
    if (_NULL_DELTA == key_bucket->_next_delta) {
      home_bucket->_first_delta = _NULL_DELTA;
    }
    else {
      home_bucket->_first_delta =
	  (home_bucket->_first_delta + key_bucket->_next_delta);
    }
  } else {
    if (_NULL_DELTA == key_bucket->_next_delta) {
      prev_key_bucket->_next_delta = _NULL_DELTA;
    }
    else {
      prev_key_bucket->_next_delta =
	  (prev_key_bucket->_next_delta + key_bucket->_next_delta);
    }
  }

  key_bucket->_next_delta = _NULL_DELTA;
  ++(_segments[segmentFunction(hash)]._timestamp);
}

// swap the free bucket within hop range if possible
// otherwise return freeBucket as NULL
template <typename IT, typename S>
bool
HopscotchSegment<IT, S>
  ::findCloserFreeBucket
    ( const IT &key
    , Segment * segment
    , Bucket **freeBucket
    , IT & freeDist
    )
{
  const IT hash (hashFunction(key));
  IT nullkey = _NULL_KEY;
  for (IT i = 0; i < _hop_range - 1; ++i)
  {
    Bucket * swapBucket (*freeBucket - _hop_range + 1 + i);
    
    // if the candidate points to a negative hash value
    if ( (*freeBucket)->_hash - _hop_range + 1 + i < 0 )
    {
      swapBucket += _size;
    }

    // in case the swap candidate now holds an empty key
    if ( swapBucket->_key == _NULL_KEY &&
         swapBucket->_key.compare_exchange_weak
	   ( nullkey, key, std::memory_order_release
	   , std::memory_order_relaxed ) )
    {
      // Here, we need to "release" the freeBucket
      (*freeBucket)->_key.store(_NULL_KEY, std::memory_order_relaxed); 
      *freeBucket = swapBucket;
      freeDist = distance(hash, (*freeBucket)->_hash);
      return true;
    }

    IT const swapDist ( distance( hashFunction(swapBucket->_key)
			        , (*freeBucket)->_hash) );

    // check if it's possible to swap our empty space with this key
    if (swapDist < _hop_range)
    {
      swapPositions
	( segment
	, freeBucket
	, swapBucket); // swap their positions

      freeDist = distance(hash, (*freeBucket)->_hash);

      return true;
    }
  }

  *freeBucket = NULL;
  
  return false;
}

template <typename IT, typename S>
void
HopscotchSegment<IT, S>
  ::swapPositions
    ( Segment * segment
    , Bucket **freeSpace
    , Bucket *swapSpace )
{
  // need to check the swapSpace for previous deltas that link to it

  IT hashValue = hashFunction(swapSpace->_key);

  // lock the segment here in case that it differs from the other lock
  std::unique_ptr<std::lock_guard<std::mutex>> lock;
  {
    if(&(_segments[segmentFunction(hashValue)])!=segment) {
      lock = std::make_unique<std::lock_guard<std::mutex>>
	  (_segments[segmentFunction(hashValue)]._mtx);
    }
  }

  // swapSpace home bucket
  Bucket *home_bucket = &(_table[hashValue]);

  // Find the previous bucket linking to swapSpace
  Bucket *prev_bucket = NULL;
  Bucket *search_bucket(home_bucket);
  IT next_delta = search_bucket->_first_delta;
  while ( (search_bucket + next_delta) != swapSpace
      &&  (search_bucket + next_delta - _size) != swapSpace
      &&  (search_bucket + next_delta + _size) != swapSpace)
  {
    advance(&search_bucket, next_delta);
    prev_bucket = search_bucket;
    next_delta = search_bucket->_next_delta;
  }

  Bucket *next_bucket = NULL;

  if(swapSpace->_next_delta != _NULL_DELTA) {
    advance(&next_bucket, next_delta);
  }
    
  
  if (NULL == prev_bucket) {
    home_bucket->_first_delta = distance(hashValue,(*freeSpace)->_hash);
  } else {
    prev_bucket->_next_delta = distance(prev_bucket->_hash,(*freeSpace)->_hash);
  }

  if(NULL == next_bucket) {
    (*freeSpace)->_next_delta = _NULL_DELTA;
  } else {
    (*freeSpace)->_next_delta = swapSpace->_next_delta
	                      - distance(swapSpace->_hash,(*freeSpace)->_hash);
  }

  // swap the element to the empty spot
  (*freeSpace)->_key.store(swapSpace->_key, std::memory_order_relaxed);
  (*freeSpace)->_value = swapSpace->_value;

  // empty out the swapped element
  {
    swapSpace->_key.store(_NULL_KEY, std::memory_order_relaxed);
    swapSpace->_value = _NULL_VALUE; 
    swapSpace->_next_delta  = _NULL_DELTA;
  }

  // update freeSpace pointer
  *freeSpace = swapSpace;
  // update time stamp in case
  if(static_cast<bool>(lock)) {
    ++(_segments[segmentFunction(hashValue)]._timestamp);
  }
}

template <typename IT, typename S>
void
HopscotchSegment<IT, S>
  ::addToNeighborhood
    ( const IT key
    , Bucket *freeBucket
    , Bucket *prevBucket
    , const S value)
{
  IT homeHashValue (hashFunction(key));

  // update key and value in the empty bucket
  freeBucket->_key.store(key, std::memory_order_relaxed);
  freeBucket->_value = value;

  // if there's a previous entry that hashes to the same value, so we need to
  // update the next deltas
  if (prevBucket != NULL)
  {
    IT const difference(distance( prevBucket->_hash
				, freeBucket->_hash));

    // if the previous bucket doesn't link to anything further, just update
    // it's next delta
    if (prevBucket->_next_delta == _NULL_DELTA)
    {
      prevBucket->_next_delta = difference;
    }
    else
    {
       // there's a value past the freeBucket
       freeBucket->_next_delta = prevBucket->_next_delta - difference;
       prevBucket->_next_delta = difference;
    }
  }
  else
  {
    // there are no previous entries that hash to this bucket, so we update the
    // first delta to the new entry
    Bucket * const start_bucket(&(_table[homeHashValue]));

    IT const difference(distance( start_bucket->_hash
    				, freeBucket->_hash));
    if (start_bucket->_first_delta == _NULL_DELTA)
    {
      start_bucket->_first_delta = difference;
    }
    else
    {
      freeBucket->_next_delta = start_bucket->_first_delta - difference;
      start_bucket->_first_delta = difference;
    }
    
  }

}

template <typename IT, typename S>
void
HopscotchSegment<IT, S>
  ::advance
    ( Bucket **current_bucket
    , const IT next_delta ) const
{
  if ( (*current_bucket) - _table.get() + next_delta >= _size )
  {
    (*current_bucket) += next_delta - _size;
  }
  else if ( (*current_bucket) - _table.get() + next_delta < 0 )
  {
    (*current_bucket) += next_delta + _size;
  }
  else
  {
    (*current_bucket) += next_delta;
  }
}

// distance between two hash values
template <typename IT, typename S>
IT
HopscotchSegment<IT, S>
  ::distance
    ( const IT & hashBeg
    , const IT & hashEnd ) const
{
  IT difference (hashEnd - hashBeg);

  // check if a loop occurs - turn this into a function later
  if (difference < 0)
  {
    difference += _size;
  }

  return difference;
}


// check if an entry is in the table
template <typename IT, typename S>
bool 
HopscotchSegment<IT, S>
    ::contains(const IT &key) const
{
  return (find(key) == -1) ? false : true;
}

// return hash value of key
template <typename IT, typename S>
IT
HopscotchSegment<IT, S>
  ::find(const IT &key) const
{
  IT hashValue = hashFunction(key);

  Bucket *current_bucket(&(_table[hashValue]));
  IT next_delta = current_bucket->_first_delta;
  Segment *segment(&(_segments[segmentFunction(hashValue)]));
  IT start_timestamp = -1;

  // iterate till the linked list is exhausted
  do {
    start_timestamp = segment->_timestamp;
    while (next_delta != _NULL_DELTA)
    {
      advance(&current_bucket, next_delta); // update current bucket to next delta
      // check if the key is contained
      if (key == current_bucket->_key)
      {
	// returns true if key is contained
	return current_bucket->_value;
      }
      // change the previous delta to the next_delta of the current bucket
      next_delta = current_bucket->_next_delta;
    }
  } while (start_timestamp != segment->_timestamp);

  return -1; // key not contained
}

template <typename IT, typename S>
std::pair<IT,bool>
HopscotchSegment<IT, S>
    ::insert(const IT &key, const S &value)
{
  IT hash(hashFunction(key));
  // check if the key is already included
  IT containedHash = find(key);

  // if already included, return the pair
  if (containedHash != -1)
  {
    //throw std::invalid_argument( "Key already included.\n" );
    std::cout << "Key already included.\n";
    _table[containedHash]._value += value;
    return std::make_pair(containedHash, false);
  }

  Segment * const segment(&(_segments[segmentFunction(hash)]));
  std::lock_guard<std::mutex> lock(segment->_mtx);

  Bucket *freeBucket(NULL);
  Bucket *prevBucket(NULL);

  IT freeDist(0);

  if(linearProbe(key, &freeBucket, &prevBucket, freeDist))
  {
    do
    {
      if(freeDist < _hop_range)
      {
	 addToNeighborhood(key, freeBucket, prevBucket, value);
	 return std::make_pair(freeBucket->_hash, true);
      }
      findCloserFreeBucket(key, segment, &freeBucket, freeDist);

    } while(freeBucket != NULL);
  }
  else
  {
    throw std::invalid_argument
      ("Table is full.");
  }


  throw std::invalid_argument
    ( "Insert failed. Need to change hoprange size.");
  //return false;
}

template <typename IT, typename S>
std::pair<IT,bool>
HopscotchSegment<IT, S>
    ::set(const IT &key, const S &value)
{
  IT hash(hashFunction(key));
  // check if the key is already included
  IT containedHash = find(key);

  // if already included, return the pair
  if (containedHash != -1)
  {
    //throw std::invalid_argument( "Key already included.\n" );
    std::cout << "Key already included.\n";
    _table[containedHash]._value = value;
    return std::make_pair(containedHash, false);
  }

  Segment * const segment(&(_segments[segmentFunction(hash)]));
  std::lock_guard<std::mutex> lock(segment->_mtx);

  Bucket *freeBucket(NULL);
  Bucket *prevBucket(NULL);

  IT freeDist(0);

  if(linearProbe(key, &freeBucket, &prevBucket, freeDist))
  {
    do
    {
      if(freeDist < _hop_range)
      {
	 addToNeighborhood(key, freeBucket, prevBucket, value);
	 return std::make_pair(freeBucket->_hash, true);
      }
      findCloserFreeBucket(key, segment, &freeBucket, freeDist);

    } while(freeBucket != NULL);
  }
  else
  {
    throw std::invalid_argument
      ( "Table is full.\n" );
  }


  throw std::invalid_argument
    ( "Insert failed. Need to change hoprange size.\n" );
  //return false;
}

// remove an element from the table
template <typename IT, typename S>
bool
HopscotchSegment<IT, S>
    ::remove(const IT &key)
{
  IT hashValue = hashFunction(key);

  Bucket *const start_bucket(&(_table[hashValue]));
  Bucket *last_bucket(NULL);
  Bucket *curr_bucket(start_bucket);
  Segment *segment(&(_segments[segmentFunction(hashValue)]));
  std::lock_guard<std::mutex> lock(segment->_mtx);

  IT next_delta(curr_bucket->_first_delta);

  while(true)
  {
    if (next_delta == _NULL_DELTA)
    {
      throw std::invalid_argument( "Key not contained.\n" );
      //return false;
    }

    advance(&curr_bucket, next_delta);
    if (key == curr_bucket->_key)
    {
      remove_key(start_bucket, curr_bucket, last_bucket, hashValue);
      return true;
    }

    last_bucket = curr_bucket;
    next_delta = curr_bucket->_next_delta;
  }

  //return false;
}

// print out the hashes and keys
template <typename IT, typename S>
void
HopscotchSegment<IT, S>
    ::displayTable()
{
  std::cout << "Hash\tKey\tTimestamp\tValue\n";
  for (IT i = 0; i < _size; ++i)
  {
    std::cout << _table[i]._hash 
              << "\t" << _table[i]._key 
              << "\t" << _table[i]._value 
              << "\n";
  }
}  

// print out the linked list's deltas
template <typename IT, typename S>
void
HopscotchSegment<IT, S>
    ::displayLists()
{
  std::cout << "Hash\tKey\tFirst Delta\tNext Delta\tHash(key)\tValue\n";
  for (IT i = 0; i < _size; ++i)
  {
    std::cout << i
	      << "\t" << _table[i]._key
	      << "\t" << _table[i]._first_delta
	      << "\t\t" << _table[i]._next_delta
	      << "\t\t" << hashFunction(_table[i]._key)
        << "\t\t" << hashFunction(_table[i]._value)
	      << "\n";
  }
}  
