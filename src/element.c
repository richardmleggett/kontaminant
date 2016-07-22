/*
 * Copyright 2009-2011 Zamin Iqbal and Mario Caccamo
 * 
 * CORTEX project contacts:  
 * 		M. Caccamo (mario.caccamo@bbsrc.ac.uk) and 
 * 		Z. Iqbal (zam@well.ox.ac.uk)
 *
 * Development team: 
 *       R. Ramirez-Gonzalez (Ricardo.Ramirez-Gonzalez@bbsrc.ac.uk)
 *       R. Leggett (richard@leggettnet.org.uk)
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
 element.c -- implements the nodes of the dBruijn graph
 */

#include <stdlib.h>
#include <global.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#ifdef THREADS
#include <pthread.h>
#endif
#include <string.h>

#include <global.h>
#include <flags.h>
#include <nucleotide.h>
#include <binary_kmer.h>
#include <element.h>

long long int visited_count = 0;

void element_assign(Element * e1, Element * e2)
{
    int i;
    
	binary_kmer_assignment_operator((*e1).kmer, (*e2).kmer);

    for (i=0; i<CONTAMINANT_FIELDS; i++) {
        e1->contaminant_flags[i] = e2->contaminant_flags[i];
    }
    
#ifdef STORE_FULL_COVERAGE
    e1->coverage[0] = e2->coverage[0];
    e1->coverage[1] = e2->coverage[1];
#endif
    
    e1->flags = e2->flags;

}

boolean element_is_key(Key key, Element e, short kmer_size)
{
	if (key == NULL) {
		printf("Do not call element_is_key wth a NULL pointer. Exiting\n");
		exit(1);
	}
    
	return binary_kmer_comparison_operator(*key, e.kmer);
}

BinaryKmer *element_get_kmer(Element * e)
{
	return &(e->kmer);
}

Key element_get_key(BinaryKmer * kmer, short kmer_size, Key preallocated_key)
{

	BinaryKmer local_rev_kmer;
	binary_kmer_initialise_to_zero(&local_rev_kmer);
    
	binary_kmer_reverse_complement(kmer, kmer_size, &local_rev_kmer);

	if (binary_kmer_less_than(local_rev_kmer, *kmer, kmer_size)) {
		binary_kmer_assignment_operator(*((BinaryKmer *)
						  preallocated_key),
						local_rev_kmer);
	} else {
		binary_kmer_assignment_operator(*((BinaryKmer *)
						  preallocated_key), *kmer);
	}

	return preallocated_key;

}

void element_initialise(Element * e, BinaryKmer * kmer, short kmer_size)
{

	BinaryKmer tmp_kmer;
    int i;

	binary_kmer_initialise_to_zero(&tmp_kmer);
	binary_kmer_assignment_operator(e->kmer, *(element_get_key(kmer, kmer_size, &tmp_kmer)));

    for (i=0; i<CONTAMINANT_FIELDS; i++) {
        e->contaminant_flags[i] = 0;
    }

    e->flags = ASSIGNED;

#ifdef STORE_FULL_COVERAGE
    e->coverage[0] = 0;
    e->coverage[1] = 0;
#endif
}

void element_set_contaminant_bit(Element* e, int id)
{
    uint32_t index = id / 32;
    uint32_t bit = id - (index * 32);
    uint32_t* field;
    //uint32_t bit = 1 << id;
    
    if (id >= MAX_CONTAMINANTS) {
        printf("Error: currently, only %d contaminants supported.\n", MAX_CONTAMINANTS);
        exit(1);
    }
    
    if (index >= CONTAMINANT_FIELDS) {
        bit += FLAG_BITS_USED;
        field = &(e->flags);
    } else {
        field = &(e->contaminant_flags[index]);
    }
    
    *field = *field | (1<<bit);
    
    //printf("id=%d     index=%d     bit=%d     field=%x\n", id, index, bit, field);
    
    //e->contaminant_flags[0] = e->contaminant_flags[0] | bit;
}

uint32_t element_get_contaminant_bit(Element* e, int id)
{
    uint32_t index = id / 32;
    uint32_t bit = id - (index * 32);
    uint32_t* field;
    //uint32_t bit = 1 << id;
    
    if (id > MAX_CONTAMINANTS) {
        printf("Error: currently, only %d contaminants supported.\n", MAX_CONTAMINANTS);
        exit(1);
    }

    if (index >= CONTAMINANT_FIELDS) {
        bit += FLAG_BITS_USED;
        field = &(e->flags);
    } else {
        field = &(e->contaminant_flags[index]);
    }
    
    return (*field) & (1<<bit);
    
    //return e->contaminant_flags[0] & bit;
}

Orientation db_node_get_orientation(BinaryKmer * k, Element * e, short kmer_size)
{
	if (binary_kmer_comparison_operator(e->kmer, *k) == true) {
		return forward;
	}

	BinaryKmer tmp_kmer;

	if (binary_kmer_comparison_operator(e->kmer,
					    *(binary_kmer_reverse_complement
					      (k, kmer_size,
					       &tmp_kmer))) == true) {
		return reverse;
	}

	printf("programming error - you have called  db_node_get_orientation with a kmer that is neither equal to the kmer in this node, nor its rev comp\n");
	char tmpseq1[kmer_size];
	char tmpseq2[kmer_size]; 
	printf("Arg 1 Kmer is %s and Arg 2 node kmer is %s\n",
	       binary_kmer_to_seq(k, kmer_size, tmpseq1),
	       binary_kmer_to_seq(&(e->kmer), kmer_size, tmpseq2));
	exit(1);

}

char binary_nucleotide_to_edge(Nucleotide base)
{
    return 1 << base;
}

Orientation opposite_orientation(Orientation o)
{
   // assert(o==forward || o == reverse); //TODO: temporary assert, to expensive to have it all the time. 
    
	return o ^ 1;

}

boolean element_check_for_flag_ALL_OFF(Element * node)
{
	return node->flags == ALL_OFF;
}

void db_node_print_binary(FILE * fp, Element* node, int kmer_size)
{
	BinaryKmer kmer;
	binary_kmer_assignment_operator(kmer, *element_get_kmer(node));
	//Edges edges = 0;
	//uint32_t coverage = 0;
	
	fwrite(&kmer, sizeof(bitfield_of_64bits), NUMBER_OF_BITFIELDS_IN_BINARY_KMER, fp);
	//fwrite(&coverage, sizeof(uint32_t), 1, fp);
	//fwrite(&edges, sizeof(Edges), 1, fp);
}

void flags_action_set_flag(Flags f, Flags * db)
{
    *db = *db | f;
}

boolean flags_check_for_any_flag(Flags f, Flags * db)
{
    boolean r = (*db & f) > 0;
    return r;
}

void element_increment_coverage(Element *node, int read)
{
#ifdef STORE_FULL_COVERAGE
    node->coverage[read]++;
#else
    if (read == 0) {
        flags_action_set_flag(COVERAGE_L, &(node->flags));
    } else if (read == 1) {
        flags_action_set_flag(COVERAGE_R, &(node->flags));
    } else {
        printf("Error - bad read number in element_increment_coverage!\n");
        exit(1);
    }
#endif
}



uint32_t element_get_coverage(Element *node, int read)
{
#ifdef STORE_FULL_COVERAGE
    return node->coverage[read];
#else
    if (flags_check_for_any_flag(read==0?COVERAGE_L:COVERAGE_R, &(node->flags))) {
        return 1;
    } else {
        return 0;
    }
#endif
}

void element_get_read_coverages(Element *node, int* a, int* b)
{
#ifdef STORE_FULL_COVERAGE
    *a = node->coverage[0];
    *b = node->coverage[1];
#else
    *a = node->flags & COVERAGE_L ? 1:0;
    *b = node->flags & COVERAGE_R ? 1:0;
#endif
}



