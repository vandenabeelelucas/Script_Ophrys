#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define MAX(x,y) ( ( x > y ) ? x : y )
#define MIN(x,y) ( ( x > y ) ? y : x )
#define LENBUFFER		10000000

/* list of input files */
char file_FASTA[200];
static char buffer[LENBUFFER], seq[LENBUFFER];

FILE *ffopen(char *File, char *Mode);

/******************************************************************************/
/* Opens a file and verifies that it opened successfully               */
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
/* Reads aligned sequence files in FASTA format         */
/******************************************************************************/
void read_file_FASTA()
{
char file_blastn[200], file_query[200], file_db[200], ident[1000], ident_max[1000];
FILE *in_FASTA, *out_ALI, *out, *in_blastn;
int nb_seq, nb_complementary, qstart, qend, sstart, send, len_query, length, l_max;
float pident;
char *ps, sseqid[500];

strcpy(file_blastn,file_FASTA);
*strchr(file_blastn,'.') = 0;
strcat(file_blastn,"-blastn.tmp");
strcpy(file_query,file_FASTA);
*strchr(file_query,'.') = 0;
strcat(file_query,"-query.tmp");
strcpy(file_db,file_FASTA);
*strchr(file_db,'.') = 0;
strcat(file_db,"-db.tmp");
*ident_max = 0;
*seq = 0;
l_max = 0;
in_FASTA = ffopen(file_FASTA,"r");
while ( fgets(buffer,LENBUFFER,in_FASTA) != NULL )
	if ( *buffer == '>' )
		{
		strcpy(ident,buffer);
		fgets(buffer,LENBUFFER,in_FASTA);
		if ( (length = strlen(buffer)) > l_max )
			{
			l_max = length;
			strcpy(ident_max,ident);
			strcpy(seq,buffer);
			}
		}
	else
		{
		fprintf(stderr,"Format error in %s: the sequence should be on a single line (use fasta2ali and ali2fasta to correct that, for instance)\n",file_FASTA);
		exit(1);
		}
fclose(in_FASTA);
out = ffopen(file_db,"w");
fprintf(out,"%s%s",ident_max,seq);
fclose(out);
sprintf (buffer,"makeblastdb -in %s -dbtype nucl -parse_seqids >/dev/null",file_db);
system (buffer);

nb_complementary = 0;
nb_seq = 0;
strcpy(buffer,file_FASTA);
*strchr(buffer,'.') = 0;
strcat(buffer,".ali");
in_FASTA = ffopen(file_FASTA,"r");
out_ALI = ffopen(buffer,"w");
fprintf(out_ALI,"#reformated by fasta2complementary\n#\n");
while ( fgets(buffer,LENBUFFER,in_FASTA) != NULL )
	if ( *buffer == '>' )
		{
		strcpy(ident,buffer);
		fgets(seq,LENBUFFER,in_FASTA);
		if ( nb_seq++ == 0 )
			fprintf(out_ALI,"%s%s",ident,seq);
		else
			{
			out = ffopen(file_query,"w");
//			fprintf(out,">query\n%s",seq);
			fprintf(out,"%s%s",ident,seq);
			fclose(out);
			sprintf(buffer,"blastn -db %s -task blastn -query %s -dust no -outfmt \"6 sseqid pident qstart qend sstart send qlen length\" -out %s",file_db,file_query,file_blastn);
			system(buffer);
//fprintf(stdout,"%s\n",buffer);
//sprintf(buffer,"cat %s %s",file_query,file_blastn);
//system(buffer);
			in_blastn = ffopen(file_blastn,"r");
			while ( fgets(buffer,LENBUFFER,in_blastn) != NULL )
				if ( *buffer != '#' )
					{
					sscanf(buffer,"%s %f %d %d %d %d %d %d", sseqid, &pident, &qstart, &qend, &sstart, &send, &len_query, &length);
					break;
					}
			fclose(in_blastn);
			if ( sstart < send )
				fprintf(out_ALI,"%s%s",ident,seq);
			else
				{
				fprintf(out_ALI,"%s",ident);
				for ( ps = seq + strlen(seq) - 2 ; ps >= seq ; --ps )
					switch(*ps)
						{
						case 'A' :
							fprintf(out_ALI,"T");
							break;
						case 'C' :
							fprintf(out_ALI,"G");
							break;
						case 'G' :
							fprintf(out_ALI,"C");
							break;
						case 'T' :
							fprintf(out_ALI,"A");
							break;
						default :
fprintf(stderr,"{%c}",*ps);
							fprintf(out_ALI,"N");
							break;
						}
				fprintf(out_ALI,"\n");
				++nb_complementary;
				}
			}
		}
fclose(in_FASTA);
fclose(out_ALI);

fprintf(stdout,"%s: %d complementary sequences among %d\n",file_FASTA,nb_complementary,nb_seq);
sprintf(buffer,"rm -f %s %s %s*",file_query,file_blastn,file_db);
//system(buffer);
}

void usage(void)
{
fprintf(stderr,"USAGE: fasta2complemantary");
fprintf(stderr,"\n\t<fasta file containing the raw data for which some sequences are on the complementary strand>\n");
exit(1);
}

/****************************************************************************/
int main(int argc, char **argv)
{
if ( argc != 2 )
	usage();

strcpy(file_FASTA,argv[1]);
read_file_FASTA();
exit(0);
}
