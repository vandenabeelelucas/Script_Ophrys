#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>

#define MAX(x,y) ( ( x > y ) ? x : y )
#define MIN(x,y) ( ( x > y ) ? y : x )
#define MAXFILEALI		35000
#define LENIDENT		80
#define LENSEQ			50000000
#define MAXSPECIES		5000
#define MAXCARSEQ		2000000000
#define LENBUFFER		50000000
#define LENMOLDOM		100
#define LENFILENAME		100

/* list of input files */
char *buffer;
int nb_file_ali;
char file_ALI[MAXFILEALI][LENFILENAME];
char molecule[LENMOLDOM], domaine[LENMOLDOM], car_inconnu[100];
int type_seq, status_gap;
int nb_spec, nb_nuc_tot, nb_nuc_var, nb_nuc_info;
char drap_fixed_type_seq;
long lenseq_tot;           /* Maximum length that sequences can have */
                           /* due to MAXCARSEQ and nb_spec */
char seq[MAXCARSEQ];
char *ident, *delcom, *retain;
int mask[MAXSPECIES], new[MAXSPECIES];

static char Mois[12][10] = { "January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"};
static char Jour[7][10] = { "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"};

void libere_memoire(void);
char read_line(FILE *in, int type, int *nb_line);
char interp_seq(int nb_line);
int cherche_delcom(void);
void date(char *version);
void date_fichier(char *nom, char *version);
FILE *ffopen(char *File, char *Mode);

/******************************************************************************/
/* Displays the fatal error message and exits the program      */
/******************************************************************************/
void erreur_fatale(char *message)
{
fprintf(stderr,"The program has encountered the following FATAL ERROR :\n\n");
fprintf(stderr,"%s\n\n\n",message);
exit(1);
}

/******************************************************************************/
/* Stores the current date in the variable version                 */
/******************************************************************************/
void date(char *version)
{
struct tm *date;
time_t l_time;

time(&l_time);
date = localtime(&l_time);
sprintf(version,"%s %d %s %d at %d hours %d",Jour[date->tm_wday],date->tm_mday,Mois[date->tm_mon],1900 + date->tm_year,date->tm_hour,date->tm_min);
}

/******************************************************************************/
/* Stores the date of use of the file 'nom' in the variable 'version'           */
/******************************************************************************/
void date_fichier(char *nom, char *version)
{
struct tm *date_f;
struct stat buf;

stat(nom,&buf);
date_f = localtime(&buf.st_atime);
sprintf(version,"%d %s %d at %d hours %d",date_f->tm_mday,Mois[date_f->tm_mon],1900 + date_f->tm_year,date_f->tm_hour,date_f->tm_min);
}

/******************************************************************************/
/* Opens a file and verifies that it opened successfully                  */
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
/* Provides the numerical code for RNA/DNA                     */
/******************************************************************************/
int code_arn(char base)
{
switch(base)
    {
    case 'A' :
        return(8);
        break;
    case 'C' :
        return(4);
        break;
    case 'G' :
        return(2);
        break;
    case 'T' :
    case 'U' :
        return(1);
        break;
    case '*' :
    case '-' :
    case '$' :
        return(0);
        break;
    case ' ' :
    case 'X' :
    case '?' :
        return(16);
        break;
    case 'N' :
        return(15);
        break;
    case 'R' :
        return(10);
        break;
    case 'Y' :
        return(5);
        break;
    case 'M' :
        return(12);
        break;
    case 'W' :
        return(9);
        break;
    case 'S' :
        return(6);
        break;
    case 'K' :
        return(3);
        break;
    case 'D' :
        return(11);
        break;
    case 'H' :
        return(13);
        break;
    case 'V' :
        return(14);
        break;
    case 'B' :
        return(7);
        break;
    default  :
        return(-1);
        break;
    }
}

/******************************************************************************/
/* Provides the numerical code for proteins                     */
/******************************************************************************/
int code_prot(char aa)
{
switch(aa)
    {
    case '*' :
    case '-' :
    case '$' :
        return(0);
        break;
    case 'A' :
        return(1);
        break;
    case 'C' :
        return(2);
        break;
    case 'D' :
        return(3);
        break;
    case 'E' :
        return(4);
        break;
    case 'F' :
        return(5);
        break;
    case 'G' :
        return(6);
        break;
    case 'H' :
        return(7);
        break;
    case 'I' :
        return(8);
        break;
    case 'K' :
        return(9);
        break;
    case 'L' :
        return(10);
        break;
    case 'M' :
        return(11);
        break;
    case 'N' :
        return(12);
        break;
    case 'P' :
        return(13);
        break;
    case 'Q' :
        return(14);
        break;
    case 'R' :
        return(15);
        break;
    case 'S' :
        return(16);
        break;
    case 'T' :
        return(17);
        break;
    case 'V' :
        return(18);
        break;
    case 'W' :
        return(19);
        break;
    case 'Y' :
        return(20);
        break;
    case ' ' :
    case '?' :
    case 'X' :
        return(21);
        break;
    default  :
        return(-1);
        break;
    }
}

/******************************************************************************/
/* Reads the list of CARAC_FILE type files                   */
/******************************************************************************/
void read_name_file_ali()
{
FILE *in;

sprintf(buffer,"ls -1 *.ali >ali2paup.tmp");
system(buffer);

nb_file_ali = 0;
in = ffopen("ali2paup.tmp","r");
while ( fgets(buffer,LENBUFFER - 1,in) != NULL )
	{
	if ( nb_file_ali >= MAXFILEALI )
		{
		fprintf(stderr,"There are too many files for this program (MAX=%d)\n",MAXFILEALI);
		exit(1);
		}
	*strchr(buffer,'\n') = 0;
	strncpy(file_ALI[nb_file_ali++],buffer,LENFILENAME);
	}
fclose(in);
system("\"rm\" ali2paup.tmp");

if ( nb_file_ali == 0 )
	{
	fprintf(stderr,"There is no file for this program\n");
	exit(1);
	}
}

/******************************************************************************/
/* Frees dynamically allocated memory                         */
/******************************************************************************/
void libere_memoire()
{
if ( ident != NULL )
	{
	free(ident);
	ident = NULL;
	}
if ( delcom != NULL )
	{
	free(delcom);
	delcom = NULL;
	}
if ( retain != NULL )
	{
	free(retain);
	retain = NULL;
	}
}

/******************************************************************************/
/* Reads a line in an .ALI file and discards comments  */
/******************************************************************************/
char read_line(FILE *in, int type, int *nb_line)
{
char *pc;

while ( fgets(buffer,LENBUFFER - 1,in) != NULL )
	{
	++*nb_line;
	if ( (pc = strchr(buffer,'\n')) == NULL )
		{
		sprintf(buffer,"Line %d is too long (MAXIMUM = %d)",*nb_line,LENSEQ - 1);
		return(-1);
		}
	if ( *buffer != '#' )
		{
		if ( (type == 0) && (*buffer != '>') )
			{
			sprintf(buffer,"Line %d is not properly formatted: need for >Identifier",*nb_line);
			return(-1);
			}
		*pc = 0;
		return(0);
		}
	}
return(-2);
}

/******************************************************************************/
/* Interprets a line containing a sequence                      */
/******************************************************************************/
char interp_seq(int nb_line)
{
char *pc;

strupr(buffer);
nb_nuc_tot = MAX(nb_nuc_tot,strlen(buffer));
if ( (nb_spec == 0) && (drap_fixed_type_seq == 0) )
	{
	if ( strpbrk(buffer,"EFILPQ") != NULL )
		type_seq = 1;
	else
		type_seq = 0;
	}

if ( type_seq == 1 )
	{
	for ( pc = buffer ; *pc != 0 ; ++pc )
		{
		if ( *pc == 'B' )
			*pc = 'X';
		if ( *pc == 'Z' )
			*pc = 'X';
		if ( strchr("ACDEFGHIKLMNPQRSTVWY*-$ X?",*pc) == NULL )
			{
			sprintf(buffer,"Character %c is not an amino acid (on line %d)",*pc,nb_line);
			return(-1);
			}
		}
	}
else
	{
	for ( pc = buffer ; *pc != 0 ; ++pc )
		{
		if ( strchr("ACGUTN*-$ RYMWSKVHDBX?",*pc) == NULL )
			{
			sprintf(buffer,"Character %c is not a base (on line %d)",*pc,nb_line);
			return(-1);
			}
		}
	}
return(0);
}

/******************************************************************************/
/* Reads aligned sequence files in ALI format         */
/******************************************************************************/
void read_file_ALI(int no_file_ali)
{
FILE *in_ali;
char *pc, *pb, *ps;
int i, j, k, nb_line, drap_erreur;

in_ali = ffopen(file_ALI[no_file_ali],"r");
fgets(buffer,LENBUFFER - 1,in_ali);
if ( *buffer != '#' )
	{
	fprintf(stderr,"First line does not begin with #\n");
	exit(1);
	}
*strchr(buffer,'\n') = 0;
strncpy(molecule,buffer + 1,LENMOLDOM - 1);
fgets(buffer,LENBUFFER - 1,in_ali);
if ( *buffer != '#' )
	{
	fprintf(stderr,"Second line does not begin with #\n");
	exit(1);
	}
*strchr(buffer,'\n') = 0;
strncpy(domaine,buffer + 1,LENMOLDOM - 1);
nb_line = 2;
for ( nb_line = 2 , nb_spec = 0 , nb_nuc_tot = 0 ; (drap_erreur = read_line(in_ali,0,&nb_line)) == 0 ; ++nb_spec )
	{
	if ( nb_spec >= MAXSPECIES )
		{
		fprintf(stderr,"Maximal number of species equals %d\n",MAXSPECIES);
		exit(1);
		}
	switch(read_line(in_ali,1,&nb_line))
		{
		case -2 :
			fprintf(stderr,"There is an identifier without a sequence at the end of the file\n");
			exit(1);
			break;
		case -1 :
			goto erreur;
			break;
		case 0 :
			if ( interp_seq(nb_line) == -1 )
				goto erreur;
			break;
		}
	}
if ( drap_erreur == -1 )
	goto erreur;
fclose(in_ali);

if ( nb_spec == 0 )
	{
	fprintf(stderr,"This file is empty\n");
	exit(1);
	}
lenseq_tot = MAXCARSEQ / (long) nb_spec;
if ( lenseq_tot < ((long) nb_nuc_tot + 1) )
	{
	fprintf(stderr,"There are too many characters, %ld (MAX = %ld)\n",(long) nb_spec * (long) (nb_nuc_tot + 1),(long) MAXCARSEQ);
	exit(1);
	}

if ( type_seq == 0 )
	strcpy(car_inconnu," ?XRYMWSKVHDBN*-");
else
	strcpy(car_inconnu," ?X*-");

if ( (ident = (char *) malloc(nb_spec * LENIDENT * sizeof(char))) == NULL )
	{
	fprintf(stderr,"Insufficient available memory (ident)\n");
	exit(-1);
	}

if ( (delcom = (char *) malloc((nb_nuc_tot + 2) * sizeof(char))) == NULL )
	{
	fprintf(stderr,"Insufficient available memory (delcom)\n");
	exit(-1);
	}

if ( (retain = (char *) malloc((nb_nuc_tot + 2) * sizeof(char))) == NULL )
	{
	fprintf(stderr,"Insufficient available memory (retain)\n");
	exit(-1);
	}

in_ali = ffopen(file_ALI[no_file_ali],"r");
for ( i = 0 ; fgets(buffer,LENBUFFER - 1,in_ali) != NULL ; )
	if ( *buffer == '>' )
		{
		if ( (pc = strstr(buffer,"#MASKED#")) == NULL )
			mask[i] = 0;
		else
			{
			mask[i] = -1;
			*pc = '\n';
			}
		if ( (pc = strstr(buffer,"#NEW#")) == NULL )
			new[i] = 0;
		else
			{
			new[i] = -1;
			*pc = '\n';
			}
		if ( strlen(buffer + 1) > LENIDENT )
			*(buffer + LENIDENT) = '\n';
		for ( pc = ident + i * LENIDENT , pb = buffer + 1 ; *pb != '\n' ; ++pc , ++pb )
			*pc = *pb;
		*pc = 0;
		}
	else
		if ( *buffer != '#' )
			{
            strupr(buffer);
			for ( ps = seq + (long) i * lenseq_tot , j = 0 , k = 0 , pc = buffer ; *pc != '\n' ; ++ps , ++pc , ++j )
				{
				if ( (*ps = *pc) == ' ' )
					++k;
				}
			for ( ; j < nb_nuc_tot ; ++j , ++k , ++ps )
				*ps = ' ';
			*ps = 0;
			++i;
			}
fclose(in_ali);
return;

erreur:
	fclose(in_ali);
	fprintf(stderr,"%s\n",buffer);
	exit(-1);
}

/******************************************************************************/
/* Searchs for common deletions specific to selected species */
/******************************************************************************/
int cherche_delcom()
{
char *ps, *pd;
int i;

for ( pd = delcom ; pd < delcom + nb_nuc_tot ; ++pd )
	*pd = 0;
*pd = 1;
/* delcom[] = -1 to indicate a shared deletion */
// delcom[] = 50 to indicate that >=50% of species are present
// delcom[] = 90 to indicate that >=90% of species are present
/* delcom[nb_nuc_tot] = 1 to indicate the end of sequences */
for ( i = 0 ; i < nb_spec ; ++i )
	for ( ps = seq + i * lenseq_tot , pd = delcom ; pd < delcom + nb_nuc_tot ; ++pd , ++ps )
		if ( strchr("* ?$-",*ps) == NULL )
			++*pd;
for ( pd = delcom ; pd < delcom + nb_nuc_tot ; ++pd )
	if ( *pd == 0 )
		*pd = -1;
	else
		if ( ((float) *pd / (float) nb_spec) >= 0.9 )
			*pd = 90;
		else
			if ( ((float) *pd / (float) nb_spec) >= 0.5 )
				*pd = 50;
			else
				*pd = 0;
//for ( pd = delcom ; pd < delcom + nb_nuc_tot ; ++pd )
//	fprintf(stdout,"%d-",*pd);
return(0);
}

/******************************************************************************/
/* Calculates indices for all sites in the sequence            */
/******************************************************************************/
void compute_stat_positions()
{
char *ps1, *ps2, *pd, *pr;
int i, j, n, nb_char, *pn, nbpres[30], nb_state1, nb_state2;

if ( type_seq == 0 )
    nb_char = 16;
else
    nb_char = 21;

nb_nuc_var = nb_nuc_info = 0;
for ( ps1 = seq , pd = delcom , pr = retain ; *ps1 != 0 ; ++ps1 , ++pd , ++pr )
    if ( *pd != -1 )
        {
/*  Fills in the table with the presence of each nucleotide                  */
        for ( pn = nbpres ; pn <= nbpres + nb_char ; ++pn )
            *pn = 0;
        if ( type_seq == 0 )
            {
            for ( ps2 = ps1 , i = 0 ; i < nb_spec ; ++i , ps2 += lenseq_tot )
                {
                if ( (n = code_arn(*ps2)) == -1 )
                    {
                    sprintf(buffer,"This character %c is unacceptable",*ps2);
                    erreur_fatale(buffer);
                    }
                ++nbpres[n];
                }
            }
        else
            {
            for ( ps2 = ps1 , i = 0 ; i < nb_spec ; ++i , ps2 += lenseq_tot )
                {
                if ( (n = code_prot(*ps2)) == -1 )
                    {
					if ( strchr("BZ",*ps2) == NULL )
						{
						sprintf(buffer,"This character %c is unacceptable",*ps2);
						erreur_fatale(buffer);
						}
					else
						{
						fprintf(stderr,"WARNING: This character %c should be replaced by X",*ps2);
						*ps2 = 'X';
						n = code_prot(*ps2);
						}
                    }
                ++nbpres[n];
                }
            }

        nb_state1 = 0;
        nb_state2 = 0;
        for ( pn = nbpres , j = 0 ; j < nb_char ; ++pn , ++j )
            {
/* Calculates the number of nucleotide types present at a site                 */
			if ( *pn >= 1 )
				{
				if ( type_seq == 0 )
					{
					if ( status_gap == 0 )
						{
						if ( (j == 0) || (j == 1) || (j == 2) || (j == 4) || (j == 8) )
							++nb_state1;
						}
					else
						{
						if ( (j == 1) || (j == 2) || (j == 4) || (j == 8) )
							++nb_state1;
						}
					}
				else
					{
					if ( status_gap == 0 )
						++nb_state1;
					else
						if ( j !=  0 )
							++nb_state1;
					}
				}
/* Calculates the number of nucleotide types present at least twice          */
			if ( *pn >= 2 )
				{
				if ( type_seq == 0 )
					{
					if ( status_gap == 0 )
						{
						if ( (j == 0) || (j == 1) || (j == 2) || (j == 4) || (j == 8) )
							++nb_state2;
						}
					else
						{
						if ( (j == 1) || (j == 2) || (j == 4) || (j == 8) )
							++nb_state2;
						}
					}
				else
					{
					if ( status_gap == 0 )
						++nb_state2;
					else
						if ( j !=  0 )
							++nb_state2;
					}
				}
/*		    if ( *pn >= 2 )
		        {
		        if ( type_seq == 0 )
		            {
		            if ( (j == 0) || (j == 1) || (j == 2) || (j == 4) || (j == 8) )
		                ++nb_state2;
		            }
		        else
		            ++nb_state2;
		        }
*/
		    }
	if ( nb_state1 > 1 )
	    {
	    ++nb_nuc_var;
	    if ( nb_state2 > 1 )
	    	{
	    	*pr = 2;
	        ++nb_nuc_info;
	    	}
	    else
	    	*pr = 1;
	    }
	else
		*pr = 0;
	}
*pr = -1;
}

/****************************************************************************/
void write_split_file_ALI(int no_file_ali)
{
char *pr, filename[LENFILENAME];
FILE *out_half1, *out_half2;
int i, j, n, n_half;

for ( pr = retain , i = 0 , n = 0 ; *pr != -1 ; ++pr , ++i )
	if ( *pr == 2 )
		{
		if ( ++n >= (nb_nuc_info / 2) )
			break;
		}
for ( j = i + 1 , ++pr ; *pr != -1 ; ++pr , ++j )
	if ( *pr == 2 )
		break;
n_half = (i + j) / 2 + 1;
fprintf(stdout,"%s has been split: first part %d (%d) positions and second part %d (%d) positions\n",file_ALI[no_file_ali],n_half,n,(int) nb_nuc_tot - n_half,nb_nuc_info - n);

strcpy(filename,file_ALI[no_file_ali]);
*strchr(filename,'.') = 0;
strcat(filename,"-half1.ali");
out_half1 = ffopen(filename,"w");
strcpy(filename,file_ALI[no_file_ali]);
*strchr(filename,'.') = 0;
strcat(filename,"-half2.ali");
out_half2 = ffopen(filename,"w");
fprintf(out_half1,"#%s\n#half1\n",file_ALI[no_file_ali]);
fprintf(out_half2,"#%s\n#half2\n",file_ALI[no_file_ali]);
for ( i = 0 ; i < nb_spec ; ++i )
	{
	*buffer = *(seq + i * lenseq_tot + n_half);
	*(seq + i * lenseq_tot + n_half) = 0;
	fprintf(out_half1,">%s\n%s\n",ident + i * LENIDENT,seq + i * lenseq_tot);
	*(seq + i * lenseq_tot + n_half) = *buffer;
	fprintf(out_half2,">%s\n%s\n",ident + i * LENIDENT,seq + i * lenseq_tot + n_half);
	}
fclose(out_half1);
fclose(out_half2);
}

/****************************************************************************/

void usage(void)
{
fprintf(stderr,"USAGE: split-gene2equal-info");
fprintf(stderr,"\n\tgap=<known or unknown> (default=unknown)");
fprintf(stderr,"\n\ttype=<protein or DNA> (optional)");
fprintf(stderr,"\n\tali=<name of file .ali> (optional)\n\n");
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

ident = NULL;
delcom = NULL;
retain = NULL;

nb_file_ali = 0;
drap_fixed_type_seq = 0;
status_gap = 1;
for( i = 1 ; i < argc ; ++i )
	{
	if ( strncmp(argv[i],"ali=",4) == 0 )
		{
		strcpy(file_ALI[0],argv[i] + 4);
		nb_file_ali = 1;
		if ( strchr(file_ALI[0],'.') != NULL )
			*strchr(file_ALI[0],'.') = 0;
		strcat(file_ALI[0],".ali");
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
				if ( strncmp(argv[i],"gap=known",9) == 0 )
					status_gap = 0;
				else
					if ( strncmp(argv[i],"gap=unknown",11) == 0 )
						status_gap = 1;
					else
						{
						fprintf(stderr,"This argument is not correctly formatted: %s\n\n",argv[i]);
						usage();
						}
	}

if ( nb_file_ali == 0 )
	read_name_file_ali();

for ( i = 0 ; i < nb_file_ali ; ++i )
	{
	read_file_ALI(i);
	cherche_delcom();
	compute_stat_positions();
	write_split_file_ALI(i);
	libere_memoire();
	}
exit(0);
}
