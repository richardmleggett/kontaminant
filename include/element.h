/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo 
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * **********************************************************************
 *
 * This file is part of CORTEX.
 *
 * CORTEX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CORTEX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CORTEX.  If not, see <http://www.gnu.org/licenses/>.
 *
 * **********************************************************************
 */

/*
 element.h defines the interface for the de Bruijn graph node. The implementation is complemented by
 a hash table that stores every node indexed by kmers (BinaryKmers).

 The element routines, ie the one required by hash_table/priority queue, are prefixed with element_
 The de Bruijn based routines are prefixed with db_node
 */

//#include <global.h>
//#include <nucleotide.h>
//#include <seq.h>
//#include <binary_kmer.h>
//#include <flags.h>

#define MAX_CONTAMINANTS 32
#define CONTAMINANT_FIELDS 1

typedef char Edges;

// To pack it, use
// typedef struct  __attribute__((packed)) {
typedef struct  __attribute__((packed)) {
 	BinaryKmer kmer;
    uint32_t contaminant_flags[CONTAMINANT_FIELDS];
    uint16_t coverage[2];
	uint16_t flags;
} Element;

typedef BinaryKmer* Key;

void element_initialise(Element * e, BinaryKmer * kmer, short kmer_size);
void element_assign(Element * e1, Element * e2);
boolean element_is_key(Key key, Element e, short kmer_size);
boolean element_check_for_flag_ALL_OFF(Element * node);
void element_set_contaminant_bit(Element* e, int id);
uint32_t element_get_contaminant_bit(Element* e, int id);
Key element_get_key(BinaryKmer * kmer, short kmer_size, Key preallocated_key);
BinaryKmer *element_get_kmer(Element * e);
void db_node_print_binary(FILE * fp, Element* node, int kmer_size);