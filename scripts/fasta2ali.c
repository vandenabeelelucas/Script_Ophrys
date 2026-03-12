#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <malloc.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>

#define MAX(x,y) ( ( x > y ) ? x : y )
#define MIN(x,y) ( ( x > y ) ? y : x )
#define MAXFILEFASTA	100000
#define LENBUFFER		50000000

/* list of presentation files */
char *buffer;
int nb_file_fasta;
char file_FASTA[MAXFILEFASTA][100], seq[LENBUFFER];
int type_seq, drap_fixed_type_seq;

FILE *ffopen(char *File, char *Mode);

/******************************************************************************/
/* Opens a file and verifies that it opened successfully                   */
/******************************************************************************/
FILE *ffopen(char *File, char *Mode)
{
FILE *in;

if ( (in = fopen(File,Mode)) == NULL )
	{
	fprintf(stderr,"Problem in opening file %s\n\n",File);
	exit(1);
	}
else
	return(in);
}

/******************************************************************************/
void strupr(char *chaine)
{
char *pc;

for ( pc = chaine ; *pc != 0 ; ++pc )
    *pc = toupper(*pc);
}

/******************************************************************************/
/* Reads the list of CARAC_FILE type files                   */
/******************************************************************************/
void read_name_file_fasta()
{
FILE *in;

sprintf(buffer,"ls -1 *.fasta >ali2paup.tmp");
system(buffer);

nb_file_fasta = 0;
in = ffopen("ali2paup.tmp","r");
while ( fgets(buffer,1000,in) != NULL )
	{
	if ( nb_file_fasta >= MAXFILEFASTA )
		{
		fprintf(stderr,"There are too many files for this program (MAX=%d)\n",MAXFILEFASTA);
		exit(1);
		}
	*strchr(buffer,'\n') = 0;
	strcpy(file_FASTA[nb_file_fasta++],buffer);
	}
fclose(in);
system("\"rm\" ali2paup.tmp");

if ( nb_file_fasta == 0 )
	{
	fprintf(stderr,"There is no file for this program\n");
	exit(1);
	}
}

/******************************************************************************/
void check_sequence(int no_file_fasta)
{
char *pc, d;

//fprintf(stdout,"{%s}\n",seq); fflush(stdout);
strupr(seq);
if ( type_seq == 1 )
	{
	for ( pc = seq ; *pc != 0 ; ++pc )
		{
		if ( *pc == 'B' )
			*pc = 'X';
		if ( *pc == 'Z' )
			*pc = 'X';
		if ( strchr("ACDEFGHIKLMNPQRSTVWY*-$ X?",*pc) == NULL )
			{
			fprintf(stderr,"Character %c is not an amino acid (%s)\n",*pc,file_FASTA[no_file_fasta]);
			*pc = 'X';
			}
		}
	}
else
	{
	for ( pc = seq ; *pc != 0 ; ++pc )
		{
		if ( strchr("ACGUTN*-$ RYMWSKVHDBX?",*pc) == NULL )
			{
			fprintf(stderr,"Character %c is not a nucleotide (%s)\n",*pc,file_FASTA[no_file_fasta]);
			*pc = 'N';
			}
		}
	}
for ( pc = seq , d = 0 ; *pc != 0 ; ++pc )
	{
	if ( *pc == '-' )
		*pc = '*';
	if ( d == 0 )
		{
		if ( *pc == '*' )
			*pc = ' ';
		else
			d = -1;
		}
	}
for ( --pc ; *pc == '*' ; --pc )
	*pc = ' ';
//fprintf(stdout,"{%s}\n",seq); fflush(stdout);
}

/******************************************************************************/
void check_identifier(int no_file_fasta)
{
char *pc;

*strchr(buffer,'\n') = 0;
if ( (pc = strchr(buffer,' ')) == NULL )
	{
	if ( strchr(buffer,'_') != NULL )
		*strchr(buffer,'_') = ' ';
	else
		strcat(buffer," sp");
	}
else
	{
	while ( strchr(pc + 1,' ') != NULL )
		*strchr(pc + 1,' ') = '_';
	}
if ( strchr(buffer,'@') == NULL )
	strcat(buffer,"@auto");
}

/******************************************************************************/
/* Reads aligned sequence files in FASTA format         */
/******************************************************************************/
void FASTA2ALI(int no_file_fasta)
{
FILE *in_fasta, *out_ali;
char *ps;

if ( drap_fixed_type_seq == 0 )
	type_seq = -1;
strcpy(buffer,file_FASTA[no_file_fasta]);
if ( strchr(buffer,'.') != NULL )
	*strchr(buffer,'.') = 0;
strcat(buffer,".ali");
fprintf(stdout,"Reformating %s into %s\n",file_FASTA[no_file_fasta],buffer);
in_fasta = ffopen(file_FASTA[no_file_fasta],"r");
out_ali = ffopen(buffer,"w");
fprintf(out_ali,"#reformatted by fasta2ali\n#\n");
*seq = 0;
while ( fgets(buffer,LENBUFFER,in_fasta) != NULL )
	{
//fprintf(stdout,"%s",buffer); fflush(stdout);
	if ( *buffer == '#' )
		fprintf(out_ali,"%s",buffer);
	else
		{
		if ( *buffer == '>' )
			{
			if ( *seq != 0 )
				{
				if ( type_seq == -1 )
					{
					if ( strpbrk(seq,"EFILPQ") != NULL )
						type_seq = 1;
					else
						type_seq = 0;
					}
				check_sequence(no_file_fasta);
				fprintf(out_ali,"%s\n",seq);
				for ( ps = seq ; *ps != 0 ; ++ps )
					*ps = 0;
//				*seq = 0;	//
				}
			check_identifier(no_file_fasta);
			fprintf(out_ali,"%s\n",buffer);
			}
		else
			{
			*strchr(buffer,'\n') = 0;
			strcat(seq,buffer);
			}
		}
	}
check_sequence(no_file_fasta);
fprintf(out_ali,"%s\n",seq);
fclose(in_fasta);
fclose(out_ali);

//sprintf(buffer,"\"mv\" %s %s",file_FASTA[no_file_fasta],file_FASTA[no_file_fasta]);
//*(buffer + strlen(buffer) - 3) = 0;
//strcat(buffer,"bak");
//system(buffer);
//sprintf(buffer,"\"mv\" file_FASTA %s",file_FASTA[no_file_fasta]);
//system(buffer);
}

/******************************************************************************/
void usage(void)
{
fprintf(stderr,"USAGE: fasta2ali");
fprintf(stderr,"fasta2ali fasta=<name of file under the fasta format (.fasta)> (optional) type=<protein or DNA> (optional)\n\t");
exit(1);
}

/****************************************************************************/
int main(int argc, char **argv)
{
int i;

if ( (buffer = (char *) malloc(LENBUFFER * sizeof(char))) == NULL )
	{
	fprintf(stderr,"Insufficient available memory\n\n");
	exit(1);
	}

nb_file_fasta = 0;
drap_fixed_type_seq = 0;
type_seq = -1;
for( i = 1 ; i < argc ; ++i )
	if ( strncmp(argv[i],"fasta=",6) == 0 )
		{
		strcpy(file_FASTA[0],argv[i] + 6);
		nb_file_fasta = 1;
		}
	else
		if ( strncmp(argv[i],"type=DNA",8) == 0 )
			{
			drap_fixed_type_seq = 1;
			type_seq = 0;
			}
		else
			if ( strncmp(argv[i],"type=protein",11) == 0 )
				{
				drap_fixed_type_seq = 1;
				type_seq = 1;
				}
			else
				usage();
if ( nb_file_fasta == 0 )
	read_name_file_fasta();

for ( i = 0 ; i < nb_file_fasta ; ++i )
	FASTA2ALI(i);
exit(0);
}
