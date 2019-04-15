/*
 * CNRSIM
 * parse_frequency.h
 * Library that generates a discrete
 * distribution of alternatives
 * alleles given the allele frequency.
 *
 * @author Riccardo Massidda
 */

#ifndef PARSE_FREQUENCY
#define PARSE_FREQUENCY

/*
 * Linear distribution.
 *
 * @param n number of alleles
 * @param p pointer to be used
 * @return p pointer to the array containing the distribution
 */
double * linear ( int n, double * p );

/*
 * Parser of dbSNP string that
 * contains allelic frequency
 *
 * @param n number of alleles
 * @param freq pointer to dbSNP string
 * @param p pointer to be used
 * @return p pointer to the array containing the distribution
 */
double * parse_db_snp_freq ( int n, char * freq, double * p );

/*
 * Parser of the allelic frequency
 * as defined in the VCF file format
 * specification
 *
 * @param n number of alleles
 * @param freq pointer to the VCF frequency array
 * @param p pointer to be used
 * @return p pointer to the array containing the distribution
 */
double * parse_af ( int n, float * freq, double * p );

#endif
