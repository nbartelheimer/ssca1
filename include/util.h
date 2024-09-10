/*
   This file is part of SSCA1.

   Copyright (C) 2008-2015, UT-Battelle, LLC.

   This product includes software produced by UT-Battelle, LLC under Contract
   No. DE-AC05-00OR22725 with the Department of Energy.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the New BSD 3-clause software license (LICENSE).

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   LICENSE for more details.

   For more information please contact the SSCA1 developers at:
   bakermb@ornl.gov
*/

#ifndef __UTIL_H
#define __UTIL_H

#include <pairwise_align.h>
#include <stdio.h>
#include <types.h>

#ifdef USE_MPI
#include <mpi.h>
extern MPI_Comm world;
extern MPI_Win window;
extern void* window_base;
extern size_t window_size;
extern void* next_window_address;
extern MPI_Request request;
#elif USE_SHMEM
#include <shmem.h>
#elif USE_GASPI
#include <GASPI.h>
#include <GASPI_Ext.h>
#include <assert.h>
#include <string.h>
extern gaspi_segment_id_t segment_id;
extern gaspi_queue_id_t queue;
extern gaspi_size_t segment_size;
extern gaspi_pointer_t segment_base;
extern gaspi_pointer_t next_segment_address;
extern gaspi_pointer_t message_buffer;
extern gaspi_offset_t loc_offset;

#define GASPI_CHECK(stmt)                                                      \
  do                                                                           \
  {                                                                            \
    const gaspi_return_t gaspi_err = (stmt);                                   \
    if(gaspi_err != GASPI_SUCCESS)                                             \
    {                                                                          \
      fprintf(stderr, "[%s:%d]: GASPI call failed with: %s \n", __FILE__,       \
              __LINE__, gaspi_error_str(gaspi_err));                                            \
      fflush(stderr);                                                                         \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
    assert(gaspi_err == GASPI_SUCCESS);                                        \
  } while(0)
#endif

#ifdef USE_MPI
#define SHORT_GET(target, source, num_elems, rank)                             \
  MPI_Get(target, num_elems, MPI_SHORT, rank, (void*)source - window_base,     \
          num_elems, MPI_SHORT, window);                                       \
  QUIET()
#elif USE_SHMEM
#define SHORT_GET(target, source, num_elems, pe)                               \
  shmem_short_get(target, source, num_elems, pe)
#elif USE_GASPI
#define SHORT_GET(target, source, num_elems, rank)                             \
  GASPI_CHECK(gaspi_read(segment_id, loc_offset, rank,       \
                         segment_id, (gaspi_pointer_t)source - segment_base,             \
                         num_elems * sizeof(short), queue, GASPI_BLOCK));      \
  QUIET(); \
  memcpy((void*)target,(void*)message_buffer,num_elems * sizeof(short));
#endif
#ifdef USE_MPI
#define SHORT_GET_NB(target, source, num_elems, rank)                          \
  MPI_Get(target, num_elems, MPI_SHORT, rank, (void*)source - window_base,     \
          num_elems, MPI_SHORT, window)
#elif USE_SHMEM
#define SHORT_GET_NB(target, source, num_elems, pe)                            \
  shmem_short_get_nbi(target, source, num_elems, pe)
#elif USE_GASPI
#define SHORT_GET_NB(target, source, num_elems, rank)                          \
  GASPI_CHECK(gaspi_read(segment_id, (gaspi_pointer_t)source - segment_base, rank,       \
                         segment_id, (gaspi_pointer_t)target - segment_base,             \
                         num_elems * sizeof(short), queue, GASPI_TEST))
#endif

#ifdef USE_MPI
#define LONG_GET(target, source, num_elems, rank)                              \
  MPI_Get(target, num_elems, MPI_LONG, rank, (void*)source - window_base,      \
          num_elems, MPI_LONG, window);                                        \
  QUIET()
#elif USE_SHMEM
#define LONG_GET(target, source, num_elems, pe)                                \
  shmem_long_get((long*)target, (long*)source, num_elems, pe)
#elif USE_GASPI
#define LONG_GET(target, source, num_elems, rank)                              \
  GASPI_CHECK(gaspi_read(segment_id, loc_offset, rank,       \
                         segment_id, (gaspi_pointer_t)source - segment_base,             \
                         num_elems * sizeof(long), queue, GASPI_BLOCK));       \
  QUIET(); \
  memcpy((void*)target,(void*)message_buffer,num_elems * sizeof(long));
#endif

#ifdef USE_MPI
#define GETMEM(target, source, length, rank)                                   \
  MPI_Get(target, length, MPI_BYTE, rank, (void*)source - window_base, length, \
          MPI_BYTE, window);                                                   \
  QUIET()
#elif USE_SHMEM
#define GETMEM(target, source, length, pe)                                     \
  shmem_getmem(target, source, length, pe)
#elif USE_GASPI
#define GETMEM(target, source, length, rank)                                   \
  GASPI_CHECK(gaspi_read(segment_id, loc_offset, rank,       \
                         segment_id, (gaspi_pointer_t)source - segment_base, length,     \
                         queue, GASPI_BLOCK));                                 \
  QUIET(); \
  memcpy((void*)target,(void*)message_buffer,length);
#endif

#ifdef USE_MPI
#define SHORT_PUT(target, source, num_elems, rank)                             \
  MPI_Put(source, num_elems, MPI_SHORT, rank, (void*)target - window_base,     \
          num_elems, MPI_SHORT, window);                                       \
  QUIET()
#elif USE_SHEM
#define SHORT_PUT(target, source, num_elems, pe)                               \
  shmem_short_put(target, source, num_elems, pe)
#elif USE_GASPI
#define SHORT_PUT(target, source, num_elems, rank)                             \
  memcpy((void*)message_buffer,(void*)source,num_elems*sizeof(short)); \
  GASPI_CHECK(gaspi_write(segment_id, loc_offset, rank,      \
                          segment_id, (gaspi_pointer_t)target - segment_base,            \
                          num_elems * sizeof(short), queue, GASPI_BLOCK));     \
  QUIET()
#endif

#ifdef USE_MPI
#define QUIET() MPI_Win_flush_all(window)
#elif USE_SHMEM
#define QUIET() shmem_quiet()
#elif USE_GASPI
#define QUIET() GASPI_CHECK(gaspi_wait(queue, GASPI_BLOCK))
#endif

#ifdef USE_MPI
#define BARRIER_ALL()                                                          \
  QUIET();                                                                     \
  MPI_Barrier(MPI_COMM_WORLD)
#elif USE_SHMEM
#define BARRIER_ALL() shmem_barrier_all()
#elif USE_GASPI
#define BARRIER_ALL()                                                          \
  QUIET();                                                                     \
  GASPI_CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
#endif

#ifdef USE_MPI
static inline int
malloc_all(size_t size, void** address)
{
  *address = next_window_address;
  next_window_address += size;
  MPI_Barrier(MPI_COMM_WORLD);
  if(next_window_address - window_base > window_size)
  {
    printf("ran out of memory!\n");
    return -1;
  }
  else
    return 0;
}
#elif USE_SHMEM
static inline int
malloc_all(size_t size, void** address)
{
  *address = shmalloc(size);
  if(*address == NULL)
    return -1;
  else
    return 0;
}
#elif USE_GASPI
static inline int
malloc_all(size_t size, void** address)
{
  *address = next_segment_address;
  next_segment_address += size;
  GASPI_CHECK(gaspi_barrier(GASPI_GROUP_ALL, GASPI_BLOCK));
  if(next_segment_address - segment_base > segment_size)
  {
    fprintf(stderr, "ran out of memory!\n");
    return -1;
  }
  return 0;
}
#endif

#ifdef USE_MPI
#define FREE_ALL(address) /* unable to free memory like this */
#elif USE_SHMEM
#define FREE_ALL(address) shfree(address)
#elif USE_GASPI
#define FREE_ALL(address)
#endif

static inline int
global_index_to_rank(const seq_t* in, const index_t codon_index)
{
  return codon_index / in->local_size;
}

static inline int
global_index_to_local_index(const seq_t* in, const index_t codon_index)
{
  return codon_index % in->local_size;
}

static inline void
fetch_from_seq(const seq_t* in, index_t const codon_index, codon_t* out)
{
  int target_ep = global_index_to_rank(in, codon_index);
  int local_index = global_index_to_local_index(in, codon_index);
  short* typed_seq = (short*)in->sequence;
  SHORT_GET((short*)out, &(typed_seq[local_index]), 1, target_ep);
}

static inline void
fetch_from_seq_nb(const seq_t* in, index_t const codon_index, codon_t* out)
{
  int target_ep = global_index_to_rank(in, codon_index);
  int local_index = global_index_to_local_index(in, codon_index);
  short* typed_seq = (short*)in->sequence;
  SHORT_GET_NB((short*)out, &(typed_seq[local_index]), 1, target_ep);
}

static inline void
write_to_seq(const seq_t* in, const index_t codon_index, codon_t data)
{
  int target_ep = global_index_to_rank(in, codon_index);
  int local_index = global_index_to_local_index(in, codon_index);
  short* typed_seq = (short*)in->sequence;
  short typed_data = (short)data;
  SHORT_PUT(&(typed_seq[local_index]), &typed_data, 1, target_ep);
}

#define WAIT_NB() QUIET()

void distribute_rng_seed(unsigned int new_seed);
void seed_rng(int adjustment);
void touch_memory(void* mem, index_t size);
index_t scrub_hyphens(good_match_t* A, seq_t* dest, seq_t* source,
                      index_t length);
void assemble_acid_chain(good_match_t* A, char* result, seq_t* chain,
                         index_t length);
void assemble_codon_chain(good_match_t* A, char* result, seq_t* chain,
                          index_t length);
score_t simple_score(good_match_t* A, seq_t* main, seq_t* match);
seq_t* alloc_global_seq(index_t seq_size);
seq_t* alloc_local_seq(index_t seq_size);
void free_global_seq(seq_t* doomed);
void free_local_seq(seq_t* doomed);

#endif

