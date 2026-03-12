#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libbcc.h>
#include <ctype.h>
#include <time.h>
#include <math.h>

#define MAX(x,y) ( ( x > y ) ? x : y )
#define MIN(x,y) ( ( x > y ) ? y : x )
#define stdprn			stdout

#define LENIDENT		80    /* maximum length of identifiers */
#define MAXSPECIES		10000
#define LENARBRE		10000000
#define LENBUFFER		LENARBRE + 1
#define LENSEQ			100000000
#define MAXCOMMENT		50    /* maximum number of comment lines */
#define LENCOMMENT		100
#define LENFILENAME		200
#define MAXCARSEQ		500000000
#define LENMOLDOM		100
#define LIMINTERNAL		0.10	// percentage of internal branches that can be longer than the re-rooting one
#define MAXCLADE		1000
#define MAXSPECPERCLADE	1000
#define MAXSUBCLADE		10000
#define MAXNODEROOT2TIP	5000

typedef struct node {        /* structure defining a node */
	float dist;
	float size_subtree;
	float diff_subtrees;
	int no_ident_arb;
	int nb_taxon_subtree;
	int nb_taxon_no_subtree;
	int div_taxon;
	int no_clade_subtree[MAXCLADE];
	int no_clade_complement_subtree[MAXCLADE];
	int nb_uniq_species_subtree;
	int no_uniq_species_subtree[MAXSPECIES];
	int no_uniq_species_complement_subtree[MAXSPECIES];
	struct node *pere;
	struct node *gauche;
	struct node *droite;
} NODE, *P_NODE;

char file_ARB[LENFILENAME], file_ALI[LENFILENAME], file_clade[LENFILENAME], file_SPLIT[LENFILENAME];
int nb_comment;
char comment[MAXCOMMENT][LENCOMMENT];
char arbre[LENARBRE];
char buffer[LENBUFFER];
int enracine;      /* is 1 if the tree is rooted */
int nb_node;
int nb_neg;
char *pc_recurs;
int nb_multif;
int nb_ident_arb;
float percent_BL_root, min_diff_subtrees;
char ident_arb[MAXSPECIES][LENIDENT];    /* identifier table */
int nb_uniq_species;
char uniq_species[MAXSPECIES][LENIDENT];    /* table of species names - without @ and without redundancy due to e.g. _NCBI */
int no_clade_uniq_species[MAXSPECIES];
int nb_outoutparalog[MAXSPECIES];    /* table counting the number of times pnode->no_uniq_species_subtree[i] is smaller than proot->no_uniq_species_subtree[i]  */
int nb_inparalog[MAXSPECIES], nb_outinparalog[MAXSPECIES];
int nbtot_outoutparalog, nbtot_outinparalog, nbtot_inparalog, nbtot_outparalog, nb_clade_outout, no_clade_outout[MAXCLADE], nb_clade_outin, no_clade_outin[MAXCLADE];

char drap_species[MAXSPECIES];
int nb_div_subtree, nb_div_no_subtree;
int max_div, nb_taxa1, nb_taxa2;
int nb_spec_subtree;
int no_spec_subtree[MAXSPECIES];
int drap_dup, drap_split;
P_NODE proot;            /* pointer to the root structure  */
P_NODE pcour;            /* pointer to the current node on the screen */
int nb_internal_branch;
float internal_length[MAXSPECIES], max_internal;
int min_nb_spec_split;

char molecule[LENMOLDOM], domaine[LENMOLDOM];
int type_spec_ali;
int nb_ident_ali;
int nb_nuc_tot;
long lenseq_tot;           /* Maximum length that sequences can have due to MAXCARSEQ and nb_ident_ali */
char seq[MAXCARSEQ];
char *ident_ali;
int mask[MAXSPECIES];
int nb_clade, nb_monophyletic_clade, nb_nonempty_clade;
char name_clade[MAXCLADE][LENIDENT];
int nb_species_clade[MAXCLADE];
int nb_uniq_species_per_clade[MAXCLADE];
int nb_uniq_species_per_clan[MAXCLADE];
int drap_presence_clan[MAXCLADE];
int no_clade_tree[MAXCLADE];
char species_clade[MAXCLADE][MAXSPECPERCLADE][LENIDENT];
int no_clade_ident_arb[MAXSPECIES];
int drap_monophyly_clade[MAXCLADE];
float maxdist_sistergroup, mindist_suspicious, min_dist_subclade;
int min_overlap;
int nb_monophyletic_subclade;
P_NODE p_subclade[MAXSUBCLADE];
P_NODE p_root2tip[MAXNODEROOT2TIP];
int list_ident_ali1[MAXSPECIES], list_ident_ali2[MAXSPECIES];
P_NODE p_LCA12[MAXCLADE], p_longest_branch[MAXCLADE];
float min_long_BL;

/* Definition of TREEPLOT functions */
void cree_dichotomie(void);
void parse(int level);
char traite_dist(char **pb, char **pa);
char traite_ident_arb(char **pb, char **pa);
char read_File(void);
P_NODE malloc_node(void);
P_NODE creer_node(char **pa, P_NODE Pere);
void verif_tree(P_NODE pnode);
P_NODE compute_tree(void);
void free_tree(P_NODE pnode);
void write_file_tax_int(P_NODE pnode, FILE *out);
void write_file_tax(void);
void cherche_taxon(P_NODE pnode);
void look_all_groupe(P_NODE pnode);
void save_tree_arb_int(P_NODE pnode,FILE *out);
void save_tree_arb(void);
void cherche_plus_droite(P_NODE pnode);
char root_left(void);

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
/* Opens a file and verifies that it opened successfully                    */
/******************************************************************************/
FILE *ffopen(char *File, char *Mode)
{
FILE *in;

if ( (in = fopen(File,Mode)) == NULL )
	{
	fprintf(stderr,"Problem in opening file %s\n",File);
	exit(1);
	}
else
	return(in);
}

/******************************************************************************/
void affich_avert(char *avert, int ligne)
{
fprintf(stderr,"%s\n",avert);
}

/******************************************************************************/
/* Frees dynamically allocated memory                        */
/******************************************************************************/
void libere_memoire()
{
if ( ident_ali != NULL )
	{
	free(ident_ali);
	ident_ali = NULL;
	}
}

/******************************************************************************/
void strupr(char *chaine)
{
char *pc;

for ( pc = chaine ; *pc != 0 ; ++pc )
    *pc = toupper(*pc);
}

/******************************************************************************/
/* Reads a line from an .ALI file and discards comments  */
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
/* Interprets a line containing a sequence                    */
/******************************************************************************/
char interp_seq(int nb_line)
{
char *pc;

strupr(buffer);
nb_nuc_tot = MAX(nb_nuc_tot,strlen(buffer));
if ( nb_ident_ali == 0 )
	{
	if ( strpbrk(buffer,"EFILPQ") != NULL )
		type_spec_ali = 1;
	else
		type_spec_ali = 0;
	}

if ( type_spec_ali == 1 )
	{
	for ( pc = buffer ; *pc != 0 ; ++pc )
		{
		if ( strchr("ACDEFGHIKLMNPQRSTVWY*-$ X?",*pc) == NULL )
			{
			sprintf(buffer,"Character %c is not an amino acid (on line %d)",*pc,nb_line);
			return(-1);
			}
		if ( strchr("*-$ X?",*pc) != NULL )
			{
			*pc = ' ';
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
		if ( strchr("*-$ X?N",*pc) != NULL )
			{
			*pc = ' ';
			}
		}
	}
return(0);
}

/******************************************************************************/
/* Reads aligned sequence files in ALI format         */
/******************************************************************************/
void read_file_ALI()
{
    FILE *in_ali;
    char *pc, *pb, *ps;
    int i, j, k, nb_line, drap_erreur;

    in_ali = ffopen(file_ALI,"r");
    fgets(buffer,LENBUFFER,in_ali);
    if ( *buffer != '#' ) {
        fprintf(stderr,"First line does not begin with #\n");
        exit(1);
    }
    *strchr(buffer,'\n') = 0;
    strncpy(molecule,buffer + 1,LENMOLDOM - 1);
    fgets(buffer,LENBUFFER,in_ali);
    if ( *buffer != '#' ) {
        fprintf(stderr,"Second line does not begin with #\n");
        exit(1);
    }
    *strchr(buffer,'\n') = 0;
    strncpy(domaine,buffer + 1,LENMOLDOM - 1);
    nb_line = 2;
    for ( nb_line = 2 , nb_ident_ali = 0 , nb_nuc_tot = 0 ;
          (drap_erreur = read_line(in_ali,0,&nb_line)) == 0 ;
          ++nb_ident_ali ){
        if ( nb_ident_ali >= MAXSPECIES ) {
            fprintf(stderr,"Maximal number of species equals %d\n",MAXSPECIES);
            exit(1);
        }
        switch(read_line(in_ali,1,&nb_line)) {
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
    if ( nb_ident_ali == 0 ) {
        fprintf(stderr,"The file, %s, is empty\n",file_ALI);
        exit(1);
    }
    lenseq_tot = MAXCARSEQ / (long) nb_ident_ali;
    if ( lenseq_tot < ((long) nb_nuc_tot + 1) ) {
        fprintf(stderr,"There are too many characters, %ld (MAX = %ld)\n",(long) nb_ident_ali * (long) (nb_nuc_tot + 1),(long) MAXCARSEQ);
        exit(1);
    }
    fclose(in_ali);


    if ( (ident_ali = (char *) malloc(nb_ident_ali * LENIDENT * sizeof(char))) == NULL ) {
        fprintf(stderr,"Insufficient available memory (ident_ali)\n");
        exit(-1);
    }

    in_ali = ffopen(file_ALI,"r");
    for ( i = 0 ; fgets(buffer,LENBUFFER,in_ali) != NULL ; )
        if ( *buffer == '>' ) {
            if ( (pc = strstr(buffer,"#MASKED#")) == NULL )
                mask[i] = 0;
            else {
                mask[i] = -1;
                *pc = '\n';
            }
            if ( strlen(buffer + 1) > LENIDENT )
                *(buffer + LENIDENT) = '\n';
            for ( pc = ident_ali + i * LENIDENT , pb = buffer + 1 ; *pb != '\n' ; ++pc , ++pb )
                *pc = *pb;
            *pc = 0;
        }
        else
            if ( *buffer != '#' ) {
                for ( ps = seq + (long) i * lenseq_tot , j = 0 , k = 0 , pc = buffer ;
                      *pc != '\n' ; ++ps , ++pc , ++j ) {
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
void read_file_CLADE()
{
FILE *in_CLADE;
char *pc, *pf, stop;

if ( *file_clade == 0 )
	return;

nb_clade = 0;
in_CLADE = ffopen(file_clade,"r");
while ( fgets(buffer,LENBUFFER - 1,in_CLADE) != NULL )
	{
//fprintf(stdout,"%s",buffer); fflush(stdout);
	*(pc = strchr(buffer,':')) = 0;
	strcpy(name_clade[nb_clade],buffer);
	nb_species_clade[nb_clade] = 0;
	stop = 0;
	for ( ; stop == 0 ; pc = pf )
		{
		pf = strpbrk(pc + 1,",\n");
		if ( *pf == '\n' )
			stop = -1;
		*pf = 0;
		strcpy(species_clade[nb_clade][nb_species_clade[nb_clade]++],pc + 1);
//fprintf(stdout,"pc + 1:%s|pf + 1:%s\n",pc + 1,pf + 1); fflush(stdout);
		}
	++nb_clade;
	}
fclose(in_CLADE);

/*
int i, j;
for ( i = 0 ; i < nb_clade ; ++i )
	{
	fprintf(stdout,"[%d]%s:",i,name_clade[i]);
	for ( j = 0 ; j < nb_species_clade[i] ; ++j )
		fprintf(stdout,",%s",species_clade[i][j]);
	fprintf(stdout,"\n"); fflush(stdout);
	}
*/
}

/******************************************************************************/
int get_no_clade_ident(char *ident)
{
int j, n;

strcpy(buffer,ident);
if ( strchr(buffer,'@') != NULL )
	*strchr(buffer,'@') = 0;
for ( n = 0 ; n < nb_clade ; ++n)
	{
	for ( j = 0 ; j < nb_species_clade[n] ; ++j )
		{
		if ( strncmp(buffer,species_clade[n][j],MAX(strlen(buffer),strlen(species_clade[n][j]))) == 0 )
			break;
		}
	if ( j != nb_species_clade[n] )
		break;
	}
/*
if ( n == nb_clade )
	{
	fprintf(stderr,"This species cannot be affected to a clade: %s\n",ident);
	exit(1);
	}
	*/
return(n);
}

/******************************************************************************/
void allocate_no_clade()
{
int i, n, nb_trash;

for ( n = 0 ; n < nb_clade ; ++n)
	{
	no_clade_tree[n] = 0;
	drap_monophyly_clade[n] = -1;
	}

for ( i = 0 , nb_trash = 0 ; i < nb_ident_arb ; ++i )
	{
	n = get_no_clade_ident(ident_arb[i]);
	if ( n == nb_clade )
		{
		fprintf(stderr,"This species cannot be affected to a clade: %s (%s)\n",ident_arb[i],file_ARB);
		strcpy(buffer,ident_arb[i]);
		if ( strchr(buffer,'@') != NULL )
			*strchr(buffer,'@') = 0;
		strcpy(species_clade[nb_clade][nb_trash++],buffer);
		}
	no_clade_ident_arb[i] = n;
	++no_clade_tree[n];
	}

if ( nb_trash != 0 )
	{
	strcpy(name_clade[nb_clade],"TRASH");
	drap_monophyly_clade[nb_clade] = -1;
	nb_species_clade[nb_clade] = nb_trash;
	++nb_clade;
	}

for ( n = 0 ; n < nb_clade ; ++n)
	{
	if ( no_clade_tree[n] > 1 )
		++nb_nonempty_clade;
//fprintf(stdout,"%s: %d\n",name_clade[n],no_clade_tree[n]);
	}
//for ( i = 0 ; i < nb_ident_arb ; ++i )
//	fprintf(stdout,"%s ---> %s\n",ident_arb[i],name_clade[no_clade_ident_arb[i]]);
}

/***********************************************************************/
void cree_dichotomie()
{
int i, niveau_parenthese;
char *pc, c1, c2;

buffer[strlen(buffer)-1] = 0;   // on enleve la virgule finale
strcat ( buffer, "):0,");
niveau_parenthese = 0;
for(i=strlen(buffer)-4; i>0; i--)
{
	if ( buffer[i] == ')' )
		niveau_parenthese++;
	if ( buffer[i] == '(' )
		niveau_parenthese--;
if (niveau_parenthese == 0 )
		break;
}
for ( pc = buffer + i + 1 , c1 = '(' ; *pc != 0 ; ++pc )
	{
	c2 = *pc;
	*pc = c1;
	c1 = c2;
	}
*pc = c1;
*(pc + 1) = 0;
nb_multif++;
}

/***********************************************************************/
void parse(int level)
{
int nb_virgules;

nb_virgules = 0;
while ( *pc_recurs != ')' )
	{
	strncat ( buffer, pc_recurs , 1);
	if ( *pc_recurs == ',' )
		nb_virgules++;
	if ( *pc_recurs == '(' )
		{
		pc_recurs++;
		parse( level+1 );
		}
	if ( *pc_recurs == 0 )
		{
		fprintf(stderr,"Error, the tree stops prematurely\n\n");
		exit(1);
		}
	if ( nb_virgules == 2 )
		{
		cree_dichotomie();
		nb_virgules = 1;
		}
	pc_recurs++;
	}
strcat(buffer,")");
}

/******************************************************************************/
// Removes PP (especially from bpcomp output)
void remove_posterior()
{
char *pa, *pb;

for ( pa = arbre , pb = buffer ; *pa != 0 ; ++pa , ++pb )
	{
	*pb = *pa;
	if ( *pa == ')' )
		while ( strchr(":;",*(pa + 1)) == NULL )
			++pa;
	}
*pb = 0;
if ( strcmp(arbre,buffer) != 0 )
	fprintf(stderr,"Posterior probabilities has been removed\n");
//fprintf(stdout,"arbre =%s\nbuffer=%s\n",arbre,buffer);
strcpy(arbre,buffer);
}

/******************************************************************************/
/* Checks whether a distance is negative                            */
/******************************************************************************/
char traite_dist(char **pb, char **pa)
{
char *p;

if ( **pb != ':' )
	{
	affich_avert("A colon is missing in the tree",10);
	return(-1);
	}
*((*pa)++) = *((*pb)++);
p = *pa;
while ( (isdigit(**pb) != 0) || (**pb == '.') || (**pb == '-') || (**pb == 'e') )
	*((*pa)++) = *((*pb)++);
--*pb;
**pa = 0;
if ( atof(p) < 0.0F )
	++nb_neg;
return(0);
}

/******************************************************************************/
/* Copies an identifier into a ident                           */
/******************************************************************************/
char traite_ident_arb(char **pb, char **pa)
{
char *pb2, *pi;

if ( (pb2 = strchr(*pb,':')) == NULL )
	{
	affich_avert("A colon following an identifier is missing in the tree",10);
	return(-1);
	}
*pb2 = 0;
strncpy(ident_arb[nb_ident_arb],*pb + 1,LENIDENT - 1);
/*
if ( (pc = strchr(ident_arb[nb_spec],'@')) != NULL )
	{
	if ( (pi = strchr(pc + 1,'@')) != NULL )
		{
		sscanf(pi + 1,"%d",&seq_lens[nb_spec]);
// lignes suivantes pour supprimer ce qu'il y a apres le 2eme @
//		for ( ; *pi != 0 ; ++pi )
//			*pi = '_';
		}
	else
		seq_lens[nb_spec] = -1;
	}
*/
for ( pi = ident_arb[nb_ident_arb] + strlen(ident_arb[nb_ident_arb]) - 1 ; *pi == '_' ; --pi )
	*pi = 0;
// Dans ce programme on ne change que le premier _ en espace
if ( (pi = strchr(ident_arb[nb_ident_arb],'_')) != NULL )
	*pi = ' ';
//while ( (pi = strchr(ident_arb[nb_ident_arb],'_')) != NULL )
//	*pi = ' ';
//if ( (pi = strstr(ident_arb[nb_ident_arb],"ZP ")) != NULL )
//	*(pi + 2) = '_';
sprintf(*pa,"%d",nb_ident_arb);
while ( **pa != 0 )
	++*pa;
*pb2 = ':';
*pb = pb2;
if ( ++nb_ident_arb >= MAXSPECIES )
	{
	sprintf(buffer,"There are too many species (MAXIMUM=%d)",MAXSPECIES);
	affich_avert(buffer,10);
	return(-1);
	}
return(traite_dist(pb,pa));
}

/******************************************************************************/
/* Data file read                                   */
/******************************************************************************/
char read_File_ARB()
{
FILE *in;
int nb_o, nb_f, nb_v;
char *pb, *pa;

in = ffopen(file_ARB,"r");
nb_comment = 0;
while(1)
	{
	if ( fgets(arbre,LENBUFFER - 1,in) == NULL )
		{
		affich_avert("This file has not enough lines",10);
		return(-1);
		}
	if ( (pa = strchr(arbre,'\n')) != NULL )
		*pa = 0;
	if ( *arbre == '#' )
		{
		if ( nb_comment < MAXCOMMENT )
			strncpy(comment[nb_comment++],arbre + 1,LENCOMMENT - 1);
		}
	else
		{
		for ( pa = arbre , pb = buffer ; *pa != 0 ; ++pa )
			if ( *pa != ' ' )
				*pb++ = *pa;
		while ( fgets(arbre,LENBUFFER - 1,in) != NULL )
			{
			for ( pa = arbre ; *pa != 0 ; ++pa )
				if ( (*pa != ' ') && (*pa != 10) )
					*pb++ = *pa;
			}
		*pb = 0;
		break;
		}
	}
fclose(in);

if ( strchr(buffer,';') == NULL )
	{
	affich_avert("There is no semicolon at the end of the tree",10);
	return(-1);
	}

//fprintf(stdout,"a={%s}\nb={%s}\n",arbre,buffer);
strcpy(arbre,buffer);
remove_posterior();

for ( nb_ident_arb = 0 , nb_o = 0 , nb_v = 0 , nb_f = 0 , nb_neg = 0 , pb = buffer , pa = arbre ; *pb != 0 ; ++pb )
	{
	switch(*pb)
		{
		case '(' :
			++nb_o;
			*pa++ = *pb;
			if ( *(pb + 1) != '(' )
				{
				if ( traite_ident_arb(&pb,&pa) == -1 )
					return(-1);
				}
			break;
		case ',' :
			++nb_v;
			*pa++ = *pb;
			if ( *(pb + 1) != '(' )
				{
				if ( traite_ident_arb(&pb,&pa) == -1 )
					return(-1);
				}
			break;
		case ')' :
			++nb_f;
			*pa++ = *pb++;
			if ( *pb != ';' )
				if ( traite_dist(&pb,&pa) == -1 )
					return(-1);
			break;
		default :
printf("\n{%s}\n\n",pb);
			affich_avert("The tree is not coherent",10);
			return(-1);
			break;
		}
	}
*pa = ';';

if ( nb_o != nb_f )
	{
	affich_avert("There are not as many ( as there are of )",10);
	return(-1);
	}
if ( nb_v == nb_o )
	enracine = 1;
else
	if ( (nb_o + 1) == nb_v )
		enracine = 0;
	else
		{
		nb_multif = 0;
		pc_recurs = arbre + 1;                  // Initializes recursion
		strcpy(buffer,"(");
		parse(0);
		sprintf(arbre,"%s;",buffer);
		if (nb_multif != 0)
			{
			sprintf(buffer,"There are %d multifurcations",nb_multif);
			affich_avert(buffer,10);
			}
		enracine = 1;
		}
if ( nb_neg != 0 )
	{
	sprintf(buffer,"There are %d negative distances",nb_neg);
	affich_avert(buffer,10);
	}
return(0);
}

/******************************************************************************/
/* Node creation                                           */
/******************************************************************************/
P_NODE malloc_node()
{
P_NODE p;
int i;

if ( (p = (P_NODE) malloc (sizeof(NODE))) == NULL )
	{
	affich_avert("Insufficient memory for creation of nodes",12);
	return(NULL);
	}
p->no_ident_arb = -2 * MAXSPECIES;
p->dist = 0.0F;
p->size_subtree = -1.0F;
p->diff_subtrees = -1.0F;
p->nb_taxon_subtree = 0;
p->nb_taxon_no_subtree = 0;
p->div_taxon = 0;
p->pere = NULL;
p->gauche = NULL;
p->droite = NULL;
for ( i = 0 ; i < MAXCLADE ; ++i )
	p->no_clade_subtree[i] = 0;
p->nb_uniq_species_subtree = 0;
for ( i = 0 ; i < MAXSPECIES ; ++i )
	{
	p->no_uniq_species_subtree[i] = 0;
	p->no_uniq_species_complement_subtree[i] = 0;
	}
return(p);
}

/******************************************************************************/
/* Node creation                                             */
/******************************************************************************/
P_NODE creer_node(char **pa, P_NODE Pere)
{
    P_NODE p;
    char *pa2, *pa3, c;

    if ( (p = malloc_node()) == NULL )
        return(NULL);

    if ( **pa == '(' ) {
/* On remplit le fils gauche */
        ++*pa;
        if ( (p->gauche = creer_node(pa,p)) == NULL )
            return(NULL);
        if ( **pa != ',' ) {
            affich_avert("There is a node where there is only one son",12);
            return(NULL);
        }

/* On remplit le fils droit */
        ++*pa;
        if ( (p->droite = creer_node(pa,p)) == NULL )
            return(NULL);
        if ( **pa != ')' ) {
            affich_avert("There is a node without )",12);
            return(NULL);
        }

/* On remplit le noeud */
        *pa += 2;
        for ( pa3 = *pa ; *pa3 != 0 ; ++pa3 )
            if ( (*pa3 == ',') || (*pa3 == ')') ) {
                c = *pa3;
                *pa3 = 0;
                break;
            }
        p->dist = (float) atof(*pa);
        p->no_ident_arb = -2 -1 * ++nb_node;
        *pa3 = c;
        *pa = pa3;
    }
    else {
        *(pa2 = strchr(*pa,':')) = 0;
        p->no_ident_arb = atoi(*pa);
        *pa2 = ':';
        for ( *(pa3 = pa2) = ':' ; *pa3 != 0 ; ++pa3 )
            if ( (*pa3 == ',') || (*pa3 == ')') ) {
                c = *pa3;
                *pa3 = 0;
                break;
            }
        p->dist = (float) atof(pa2 + 1);
        *pa3 = c;
        *pa = pa3;
        p->gauche = NULL;
        p->droite = NULL;
    }
    p->pere = Pere;
    return(p);
}

/******************************************************************************/
/* Tree construction                                        */
/******************************************************************************/
P_NODE compute_tree()
{
    P_NODE p1, p2, p3, p4;
    char *pa;

	nb_node = 0;
    percent_BL_root = 1.0;
    if ( (proot = malloc_node()) == NULL )
        return(NULL);
    if ( enracine == 0 ) {
        if ( (p4 = malloc_node()) == NULL )
            return(NULL);
    }
    else
        p4 = proot;

    pa = arbre + 1;
    if ( (p1 = creer_node(&pa,p4)) == NULL )
        return(NULL);
    if ( *pa != ',' ) {
        affich_avert("The first node has not a sister-group",12);
        return(NULL);
    }
    ++pa;
    if ( (p2 = creer_node(&pa,p4)) == NULL )
        return(NULL);

/* The third subtree is created only if the tree is not rooted */
    if ( enracine == 0 ) {
        if ( *pa != ',' ) {
            affich_avert("The second node has not a sister-group",12);
            return(NULL);
        }
        ++pa;
        if ( (p3 = creer_node(&pa,proot)) == NULL )
            return(NULL);
    }

    if ( *pa != ')' ) {
        affich_avert("The final parenthesis is missing",12);
        return(NULL);
    }

    if ( enracine == 0 ) {
        p4->pere = proot;
        p4->droite = p1;
        p4->gauche = p2;
        p4->no_ident_arb = -2;

        proot->pere = NULL;
        proot->gauche = p4;
        proot->droite = p3;
        proot->no_ident_arb = -1;

        p4->dist = (p3->dist) * percent_BL_root;
        p3->dist = p3->dist * (1.0 - percent_BL_root);
    }
    else {
        proot->pere = NULL;
        proot->droite = p1;
        proot->gauche = p2;
    }
    return(proot);
}

/******************************************************************************/
/* Removes a tree by freeing allocated memory                        */
/******************************************************************************/
void free_tree(P_NODE pnode)
{
    if ( (pnode->gauche != NULL) && (pnode->droite != NULL) ) {
        free_tree(pnode->gauche);
        free_tree(pnode->droite);
    }
    free(pnode);
}

/******************************************************************************/
/* Calculates the size of each subtree                        */
/******************************************************************************/
void compute_size_subtree(P_NODE pnode)
{
if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	compute_size_subtree(pnode->gauche);
	compute_size_subtree(pnode->droite);
	pnode->size_subtree = pnode->dist + (pnode->gauche)->size_subtree + (pnode->droite)->size_subtree;
	}
else
	pnode->size_subtree = pnode->dist;
}

/******************************************************************************/
/* Function normalizing the size of the tree to 100                           */
/******************************************************************************/
void normalise100_tree(P_NODE pnode, float size_tree)
{
if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	normalise100_tree(pnode->gauche,size_tree);
	normalise100_tree(pnode->droite,size_tree);
	pnode->size_subtree *= (100.0 / size_tree);
	pnode->dist *= (100.0 / size_tree);
	}
else
	{
	pnode->size_subtree *= (100.0 / size_tree);
	pnode->dist *= (100.0 / size_tree);
	}
}

/******************************************************************************/
/* Function affecting the taxonomic distribution of each node                 */
/******************************************************************************/
void compute_taxonomy_subtree(P_NODE pnode)
{
int i;

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	compute_taxonomy_subtree(pnode->gauche);
	compute_taxonomy_subtree(pnode->droite);
	for ( i = 0 ; i < nb_clade ; ++i )
		pnode->no_clade_subtree[i] = (pnode->gauche)->no_clade_subtree[i] + (pnode->droite)->no_clade_subtree[i];
	for ( i = 0 , pnode->nb_uniq_species_subtree = 0 ; i < nb_uniq_species ; ++i )
		if ( (pnode->no_uniq_species_subtree[i] = (pnode->gauche)->no_uniq_species_subtree[i] + (pnode->droite)->no_uniq_species_subtree[i]) > 0 )
			++pnode->nb_uniq_species_subtree;
	}
else
	{
	++pnode->no_clade_subtree[no_clade_ident_arb[pnode->no_ident_arb]];
	for ( i = 0 ; i < nb_uniq_species ; ++i )
		if ( strncmp(ident_arb[pnode->no_ident_arb],uniq_species[i],strlen(uniq_species[i])) == 0 )
			break;
	if ( i == nb_uniq_species )
		fprintf(stderr,"MAJOR PROBLEM with %s\n",ident_arb[pnode->no_ident_arb]);
	else
		++(pnode->no_uniq_species_subtree[i]);
	pnode->nb_uniq_species_subtree = 1;
	}
}

/******************************************************************************/
/* Function affecting the taxonomic distribution of each node                 */
/******************************************************************************/
void erase_taxonomy_subtree(P_NODE pnode)
{
int i;

for ( i = 0 ; i < nb_clade ; ++i )
	pnode->no_clade_subtree[i] = 0;
for ( i = 0 ; i < nb_uniq_species ; ++i )
	pnode->no_uniq_species_subtree[i] = 0;
if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	erase_taxonomy_subtree(pnode->gauche);
	erase_taxonomy_subtree(pnode->droite);
	}
}

/******************************************************************************/
/* Function Removing a tree by freeing allocated memory                        */
/******************************************************************************/
void verif_tree(P_NODE pnode)
{
int i;

fprintf(stdout,"\nNoeud=%3d Dist=%.2f Size=%.2f ",pnode->no_ident_arb,pnode->dist,pnode->size_subtree); fflush(stdout);
for ( i = 0 ; i < nb_uniq_species ; ++i )
	fprintf(stdout,"%d-",pnode->no_uniq_species_subtree[i]);
if ( pnode->gauche == NULL )
	{
	fprintf(stdout," Species=%s",ident_arb[pnode->no_ident_arb]); fflush(stdout);
	}

if ( pnode->pere != NULL )
	{
	if ( pnode->gauche != NULL )
		{
		fprintf(stdout," Pere=%3d Gauche=%3d Droite=%3d",(pnode->pere)->no_ident_arb,(pnode->gauche)->no_ident_arb,(pnode->droite)->no_ident_arb); fflush(stdout);
		}
	else
		{
		fprintf(stdout," Pere=%3d",(pnode->pere)->no_ident_arb); fflush(stdout);
		}
	}
else
	{
	if ( pnode->gauche != NULL )
		{
		fprintf(stdout," Gauche=%3d Droite=%3d",(pnode->gauche)->no_ident_arb,(pnode->droite)->no_ident_arb); fflush(stdout);
		}
	else
		{
		fflush(stdout);
		}
	}

if ( pnode->gauche != NULL )
	{
	verif_tree(pnode->gauche);
	verif_tree(pnode->droite);
	}
}

/******************************************************************************/
/* Recursive function to save the tree                            */
/******************************************************************************/
void save_tree_arb_int(P_NODE pnode,FILE *out)
{
char *pi, *pb;
int i;

if ( (pnode->droite == NULL) && (pnode->gauche == NULL) )
	{
	for ( pi = ident_arb[pnode->no_ident_arb] , pb = buffer , i = 0 ; i < LENIDENT - 1 ; ++pi , ++pb , ++i )
		if ( *pi == ' ' )
			*pb = '_';
		else
			*pb = *pi;
	if ( strlen(buffer) >= LENIDENT )
		buffer[LENIDENT - 1] = 0;
	else
		{
		for ( pb = buffer + strlen(buffer) ; pb < buffer + LENIDENT - 1 ; ++pb )
			*pb = '_';
		*pb = 0;
		}
	fprintf(out,"%s",buffer);
	}
else
	{
	fprintf(out,"(");
	save_tree_arb_int(pnode->gauche,out);
	fprintf(out,",");
	save_tree_arb_int(pnode->droite,out);
	fprintf(out,")");
	}
fprintf(out,":%f",pnode->dist);
}

/******************************************************************************/
/* Function saving the trees                                                  */
/******************************************************************************/
void save_tree_arb()
{
FILE *out;
int i;

out = ffopen(file_ARB,"w");
for ( i = 0 ; i < nb_comment ; ++i )
	fprintf(out,"#%s\n",comment[i]);
fprintf(out,"(");
save_tree_arb_int(proot->droite,out);
fprintf(out,",");
save_tree_arb_int(proot->gauche,out);
fprintf(out,");\n");
fclose(out);
}

/******************************************************************************/
/* Roots the tree at the left of the active node                               */
/******************************************************************************/
char root_left()
{
    P_NODE p0, p1, p2, pr;
    float f, ff;

    if ( (pcour != proot) && (pcour->pere != proot) ) {
        if ( (pr = malloc_node()) == NULL ) {
            free_tree(proot);
            proot = NULL;
            return(-1);
        }
        pr->pere = NULL;
        pr->gauche = pcour;
        pr->droite = pcour->pere;

/* Adjust the 2 child distances of proot */
        f = (proot->gauche)->dist + (proot->droite)->dist;
        (proot->gauche)->dist = f;
        (proot->droite)->dist = f;

        for ( p0 = pcour , p1 = pcour->pere , p2 = (pcour->pere)->pere , f = p1->dist ; p1 != proot ; ) {
//printf("\np0=%d p1=%d p2=%d\n",p0->no_ident_arb,p1->no_ident_arb,p2->no_ident_arb);
//printf("d0=%.2f d1=%.2f d2=%.2f\n",p0->dist,p1->dist,p2->dist);
            if ( p2 != NULL ) {
                ff = p2->dist;
                p2->dist = f;
                f = ff;
            }
//printf("d0=%.2f d1=%.2f d2=%.2f\n",p0->dist,p1->dist,p2->dist); getch();
            if ( p1->gauche == p0 )
                p1->gauche = p1->pere;
            else
                p1->droite = p1->pere;
            p1->pere = p0;
            if ( p2->gauche == p1 )
                p2->gauche = p1;
            else
                p2->droite = p1;
            p0 = p1;
            p1 = p2;
            p2 = p2->pere;
        }

        if ( p1->gauche == p0 ) {
            if ( p0->gauche == p1 ) {
                p0->gauche = p1->droite;
                (p0->gauche)->pere = p0;
            }
            else {
                p0->droite = p1->droite;
                (p0->droite)->pere = p0;
            }
        }
        else {
            if ( p0->gauche == p1 ) {
                p0->gauche = p1->gauche;
                (p0->gauche)->pere = p0;
            }
            else {
                p0->droite = p1->gauche;
                (p0->droite)->pere = p0;
            }
        }
        (pr->droite)->pere = pr;
        pcour->pere = pr;
        f = pcour->dist;
        (pr->gauche)->dist = f * (1.0 - percent_BL_root);
        (pr->droite)->dist = f * percent_BL_root;
//verif_tree(pr);

        free(proot);
        proot = pr;
        return(0);
    }
    else {
        return(0);
    }
}

/******************************************************************************/
void lookfor_internal_distance(P_NODE pnode, int *n)
{
if ( (pnode->droite != NULL) && (pnode->gauche != NULL) )
	{
	internal_length[(*n)++] = pnode->dist;
	lookfor_internal_distance(pnode->droite,n);
	lookfor_internal_distance(pnode->gauche,n);
	}
}

/******************************************************************************/
int comparer_float(float *arg1, float *arg2)
{
if ( *arg1 > *arg2 )
	return(-1);
else
	return(1);
}

/******************************************************************************/
void lookfor_node(P_NODE pnode, int i)
{
//fprintf(stdout,"Noeud=%3d Dist=%.2f Size=%.2f pcour=%d\n",pnode->no_ident_arb,pnode->dist,pnode->size_subtree,(int) pcour); fflush(stdout);
if ( pcour != NULL )
	return;
if ( (pnode->droite != NULL) && (pnode->gauche != NULL) )
	{
	lookfor_node(pnode->droite,i);
	lookfor_node(pnode->gauche,i);
	}
if ( pnode->no_ident_arb == i )
	pcour = pnode;
}

/******************************************************************************/
/* Roots at the left of the node of the species with no_ident 0    */
/******************************************************************************/
void root_left_species0()
{
pcour = NULL;
lookfor_node(proot,0);
//fprintf(stdout,"\nNoeud=%3d Dist=%.2f Size=%.2f",pcour->no_ident_arb,pcour->dist,pcour->size_subtree); fflush(stdout);
root_left();
}

/******************************************************************************/
void lookfor_suspicious_sistergroup(P_NODE pnode)
{
int i, n;

if ( pnode->droite == NULL )
	return;
//fprintf(stdout,"(%d)\n",pnode->no_ident_arb); fflush(stdout);

// case where the two sons are leaves
if (((pnode->droite)->gauche == NULL) && ((pnode->gauche)->gauche == NULL))
	{
	if (((pnode->gauche)->dist + (pnode->droite)->dist) <= maxdist_sistergroup )
		{
		strcpy(buffer,ident_arb[(pnode->gauche)->no_ident_arb]);
		if ( strchr(buffer,'@') != NULL )
			*strchr(buffer,'@') = 0;
		strcpy(buffer + 1000,ident_arb[(pnode->droite)->no_ident_arb]);
		if ( strchr(buffer + 1000,'@') != NULL )
			*strchr(buffer + 1000,'@') = 0;
		if ( strncmp(buffer,buffer + 1000,MIN(strlen(buffer),strlen(buffer + 1000))) != 0 )
			{
			if ( no_clade_ident_arb[(pnode->gauche)->no_ident_arb] == no_clade_ident_arb[(pnode->droite)->no_ident_arb] )
				{
				if ( strncmp(buffer,buffer + 1000,MIN(strlen(buffer),strlen(buffer + 1000))) <= 0 )
					fprintf(stdout,"Suspicious sister-group for the same clade for %s (dist=%.3f): %s/%s (%s,%s)\n",file_ARB,(pnode->gauche)->dist + (pnode->droite)->dist,buffer,buffer+1000,ident_arb[(pnode->gauche)->no_ident_arb],ident_arb[(pnode->droite)->no_ident_arb]);
				else
					fprintf(stdout,"Suspicious sister-group for the same clade for %s (dist=%.3f): %s/%s (%s,%s)\n",file_ARB,(pnode->gauche)->dist + (pnode->droite)->dist,buffer+1000,buffer,ident_arb[(pnode->gauche)->no_ident_arb],ident_arb[(pnode->droite)->no_ident_arb]);
				}
			else
				{
				if ( strncmp(buffer,buffer + 1000,MIN(strlen(buffer),strlen(buffer + 1000))) <= 0 )
					fprintf(stdout,"SUSPICIOUS sister-group for different clades for %s (dist=%.3f): %s/%s (%s,%s)\n",file_ARB,(pnode->gauche)->dist + (pnode->droite)->dist,buffer,buffer+1000,ident_arb[(pnode->gauche)->no_ident_arb],ident_arb[(pnode->droite)->no_ident_arb]);
				else
					fprintf(stdout,"SUSPICIOUS sister-group for different clades for %s (dist=%.3f): %s/%s (%s,%s)\n",file_ARB,(pnode->gauche)->dist + (pnode->droite)->dist,buffer+1000,buffer,ident_arb[(pnode->gauche)->no_ident_arb],ident_arb[(pnode->droite)->no_ident_arb]);
//++pnode->no_clade_subtree[no_clade_ident_arb[pnode->no_ident_arb]];
				n = no_clade_ident_arb[(pnode->gauche)->no_ident_arb];
				for ( i = 0 ; i < nb_clade ; ++i )
					if ( (pnode->pere)->no_clade_subtree[i] > 1 )
						break;
//printf("%s: %d\n",name_clade[i],(pnode->pere)->no_clade_subtree[i]);
				if ( i == n )
					fprintf(stdout,"remove_seq \"%s\" %s\t#included within %s\n",ident_arb[(pnode->droite)->no_ident_arb],file_ALI,name_clade[i]);
				else
					fprintf(stdout,"remove_seq \"%s\" %s\t#included within %s\n",ident_arb[(pnode->gauche)->no_ident_arb],file_ALI,name_clade[i]);
				}
//				fprintf(stdout,"SUSPICIOUS sister-group for different clades for %s (dist=%.3f): (%s,%s)\n",file_ARB,(pnode->gauche)->dist + (pnode->droite)->dist,ident_arb[(pnode->gauche)->no_ident_arb],ident_arb[(pnode->droite)->no_ident_arb]);
			}
		}
	}
else
	{
	if ( pnode->droite != NULL )
		lookfor_suspicious_sistergroup(pnode->droite);
	if ( pnode->gauche != NULL )
		lookfor_suspicious_sistergroup(pnode->gauche);
	}
}

/******************************************************************************/
void affich_suspicious_species(P_NODE pnode, int no_biggest_clade)
{
if ( pnode->droite == NULL )
	{
	if ( no_clade_ident_arb[pnode->no_ident_arb] != no_biggest_clade )
		fprintf(stdout,"%s: %s incorrectly located within %s\n",file_ARB,ident_arb[pnode->no_ident_arb],name_clade[no_biggest_clade]);
	}
else
	{
	affich_suspicious_species(pnode->droite,no_biggest_clade);
	affich_suspicious_species(pnode->gauche,no_biggest_clade);
	}
}

/******************************************************************************/
void lookfor_clade_suspicious_taxonomy(P_NODE pnode)
{
int i, nb_nonzero, nb_almost_complete_clade, no_biggest_clade;

/*
fprintf(stdout,"NODE%d:",pnode->no_ident_arb);
for ( i = 0 ; i < nb_clade ; ++i )
	{
fprintf(stdout," %d",pnode->no_clade_subtree[i]);
	}
fprintf(stdout,"\n");
*/
if ( pnode->droite == NULL )
	return;

for ( i = 0 , nb_nonzero = 0 , nb_almost_complete_clade = 0 , no_biggest_clade = 0 ; i < nb_clade ; ++i )
	{
	if ( pnode->no_clade_subtree[i] != 0 )
		++nb_nonzero;
	if ( (pnode->no_clade_subtree[i] == no_clade_tree[i]) && (pnode->no_clade_subtree[i] > 1) )
		++nb_almost_complete_clade;
	if ( pnode->no_clade_subtree[i] > pnode->no_clade_subtree[no_biggest_clade] )
		no_biggest_clade = i;
	}

if ( nb_nonzero == 1 )
	{
	for ( i = 0 ; i < nb_clade ; ++i )
		{
//fprintf(stdout,"%s: clade %s \n",file_ARB,name_clade[i]);
		if ( (pnode->no_clade_subtree[i] == no_clade_tree[i]) && (pnode->no_clade_subtree[i] > 1) )
			{
			fprintf(stdout,"%s: clade %s monophyletic\n",file_ARB,name_clade[i]);
			++nb_monophyletic_clade;
			drap_monophyly_clade[i] = 0;
			break;
			}
		}
	}
else
	{
	if ( nb_almost_complete_clade == 1 )
		{
//fprintf(stdout,"%s: %d %d",name_clade[no_biggest_clade],(pnode->gauche)->no_clade_subtree[no_biggest_clade],(pnode->droite)->no_clade_subtree[no_biggest_clade]);
// checking that the smallest clades are not fully embedded into this node
	for ( i = 0 ; i < nb_clade ; ++i )
		if ( (i != no_biggest_clade) && (pnode->no_clade_subtree[i] != 0) )
			if ( pnode->no_clade_subtree[i] != no_clade_tree[i] )
				break;

// checking that the biggest clade is simply one of the son
		if ( (pnode->no_clade_subtree[no_biggest_clade] != (pnode->gauche)->no_clade_subtree[no_biggest_clade]) && (pnode->no_clade_subtree[no_biggest_clade] != (pnode->droite)->no_clade_subtree[no_biggest_clade]) )
			if ( (pnode->dist >= mindist_suspicious) && (i != nb_clade) )
				{
				affich_suspicious_species(pnode,no_biggest_clade);
for ( i = 0 ; i < nb_clade ; ++i )
	fprintf(stdout,"%2d ",pnode->no_clade_subtree[i]);
fprintf(stdout,"\n");
				}
		}
	}
lookfor_clade_suspicious_taxonomy(pnode->gauche);
lookfor_clade_suspicious_taxonomy(pnode->droite);
}

/******************************************************************************/
void lookfor_monophyletic_subclade(P_NODE pnode, int i_clade)
{
int i;

for ( i = 0 ; i < nb_clade ; ++i )
	if ( i != i_clade )
		{
		if ( pnode->no_clade_subtree[i] != 0 )
			break;
		}
if ( i == nb_clade )
	{
//for ( i = 0 ; i < nb_clade ; ++i )
//	fprintf(stdout,"%2d ",pnode->no_clade_subtree[i]);
//fprintf(stdout,"\n");
	p_subclade[nb_monophyletic_subclade] = pnode;
	if ( ++nb_monophyletic_subclade >= MAXSUBCLADE )
		erreur_fatale("Too many subclades\n");
	return;
	}

if ( (pnode->droite != NULL) && (pnode->gauche != NULL) )
	{
	lookfor_monophyletic_subclade(pnode->droite,i_clade);
	lookfor_monophyletic_subclade(pnode->gauche,i_clade);
	}
}

/******************************************************************************/
void compute_dist_of_two_nodes(P_NODE p_subclade1, P_NODE p_subclade2, float *dist12, int i_clade)
{
P_NODE pnode;
int i, nb_node1R;
float max_d;

// filling the array with all the the ancestral nodes to p_subclade1
for ( pnode = p_subclade1 , nb_node1R = 0 ; pnode->pere != NULL ; pnode = pnode->pere )
	p_root2tip[nb_node1R++] = pnode;
for ( pnode = p_subclade2 ; pnode->pere != NULL ; pnode = pnode->pere )
	{
	for ( i = 0 ; i < nb_node1R ; ++i )
		if ( p_root2tip[i] == pnode )
			break;
	if ( i != nb_node1R )
		break;
	}
p_LCA12[i_clade] = pnode;
p_longest_branch[i_clade] = p_subclade1;
max_d = (p_longest_branch[i_clade])->dist;
//fprintf(stdout,"p_LCA12=%d\n",p_LCA12->no_ident_arb);
*dist12 = 0;
for ( pnode = p_subclade1 ; pnode != p_LCA12[i_clade] ; pnode = pnode->pere )
	{
	*dist12 += pnode->dist;
	if ( pnode->dist > max_d )
		{
		max_d = pnode->dist;
		p_longest_branch[i_clade] = pnode;
		}
//fprintf(stdout,"(1)%d: %f\n",pnode->no_ident_arb,pnode->dist);
	}
for ( pnode = p_subclade2 ; pnode != p_LCA12[i_clade] ; pnode = pnode->pere )
	{
	*dist12 += pnode->dist;
	if ( pnode->dist > max_d )
		{
		max_d = pnode->dist;
		p_longest_branch[i_clade] = pnode;
		}
//fprintf(stdout,"(2)%d: %f\n",pnode->no_ident_arb,pnode->dist);
	}
}

/******************************************************************************/
void lookfor_ident_ali(P_NODE pnode, int *nb_id, int *list_ident_ali)
{
int i;

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	lookfor_ident_ali(pnode->gauche,nb_id,list_ident_ali);
	lookfor_ident_ali(pnode->droite,nb_id,list_ident_ali);
	}
else
	{
	for ( i = 0 ; i < nb_ident_ali ; ++i )
		if ( strncmp(ident_arb[pnode->no_ident_arb],ident_ali + i * LENIDENT,strlen(ident_arb[pnode->no_ident_arb])) ==  0 )
			break;
	if ( i == nb_ident_ali )
		{
		fprintf(stderr,"This identifier (%s) is missing in the alignment file (%s)\n",ident_arb[pnode->no_ident_arb],file_ALI);
		exit(1);
		}
	else
		*(list_ident_ali + (*nb_id)++) = i;
	}
}

/******************************************************************************/
// function looking for the possible absence of overlap between the sequences of two sub-clades
int lookfor_overlap_subclades(P_NODE p_subclade1, P_NODE p_subclade2)
{
int i, nb_sp1, nb_sp2, nb_aa1, nb_aa2, nb_overlap;
char *ps;

nb_sp1 = nb_sp2 = 0;
lookfor_ident_ali(p_subclade1,&nb_sp1,list_ident_ali1);
lookfor_ident_ali(p_subclade2,&nb_sp2,list_ident_ali2);
//fprintf(stdout,"CLADE1 (%d sp):",nb_sp1);
//for ( i = 0 ; i < nb_sp1 ; ++i )
//	fprintf(stdout," %s",ident_ali + list_ident_ali1[i] * LENIDENT);
//fprintf(stdout,"\n");
//fprintf(stdout,"CLADE2 (%d sp):",nb_sp2);
//for ( i = 0 ; i < nb_sp2 ; ++i )
//	fprintf(stdout," %s",ident_ali + list_ident_ali2[i] * LENIDENT);
//fprintf(stdout,"\n");

for ( ps = seq , nb_overlap = 0; *ps != 0 ; ++ps )
	{
	nb_aa1 = 0;
	for ( i = 0 ; i < nb_sp1 ; ++i )
		if ( strchr("*-$ X?",*(ps + list_ident_ali1[i] * lenseq_tot)) == NULL )
			++nb_aa1;
	nb_aa2 = 0;
	for ( i = 0 ; i < nb_sp2 ; ++i )
		if ( strchr("*-$ X?",*(ps + list_ident_ali2[i] * lenseq_tot)) == NULL )
			++nb_aa2;
	if ( (nb_aa1 > 0) && (nb_aa2 > 0) )
		++nb_overlap;
	}
//fprintf(stdout,"----> %d overlapping aa\n",nb_overlap);
if ( nb_overlap > min_overlap )
	return(0);
else
	return(-1);
}

/******************************************************************************/
/* Finds the node for midpoint rooting                        */
/******************************************************************************/
void lookfor_midlength(P_NODE pnode)
{
float f, f1, f2, f3;

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	lookfor_midlength(pnode->gauche);
	lookfor_midlength(pnode->droite);
	}
f1 = fabs(proot->size_subtree - 2 * pnode->size_subtree + 2 *pnode->dist);
f2 = fabs(proot->size_subtree - 2 * pnode->size_subtree);
f = MIN(f1,f2);
f3 = fabs(proot->size_subtree - 2 * pnode->size_subtree + pnode->dist);
f = MIN(f,f3);
pnode->diff_subtrees = f;
if ( f < min_diff_subtrees )
	{
	min_diff_subtrees = f;
	pcour = pnode;
//fprintf(stdout,"Node %d: MADDIV=%d #taxa1=%d #taxa2=%d (max=%d) (diff_size= %f/%f =%f)\n",pnode->no_ident_arb,max_div,nb_taxa1,nb_taxa2,nb_uniq_species,pnode->diff_subtrees,min_diff_subtrees,100.0*pnode->diff_subtrees/min_diff_subtrees);
	}
//fprintf(stdout,"Node=%d (b=%f ecart1-2-3=%f %f %f (min=%f)\n",pnode->no_ident_arb,pnode->dist,f1,f2,f3,min_diff_subtrees);
}

/******************************************************************************/
/* Function generating the list of all species in a subtree           */
/****************************************************************************/
void lookfor_species_subtree(P_NODE pnode)
{
if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	lookfor_species_subtree(pnode->gauche);
	lookfor_species_subtree(pnode->droite);
	}
else
	no_spec_subtree[nb_spec_subtree++] = pnode->no_ident_arb;
}

/******************************************************************************/
/* Function searching in a subtree if a species is present in the no_spec_subtree list */
/****************************************************************************/
void lookfor_dup_species(P_NODE pnode)
{
int i, l;

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	lookfor_dup_species(pnode->gauche);
	lookfor_dup_species(pnode->droite);
	}
else
	{
	strcpy(buffer,ident_arb[pnode->no_ident_arb]);
	if ( strchr(buffer,'@') != NULL )
		*strchr(buffer,'@') = 0;
	l = strlen(buffer);
	for ( i = 0 ; i < nb_spec_subtree ; ++i )
		if ( strncmp(buffer,ident_arb[no_spec_subtree[i]],l) == 0 )
			{
			drap_dup = 0;
			return;
			}
	}
}

/****************************************************************************/
void generate_liste_uniq_species(P_NODE pnode)
{
int i, l;

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	generate_liste_uniq_species(pnode->gauche);
	generate_liste_uniq_species(pnode->droite);
	}
else
	{
	strcpy(buffer,ident_arb[pnode->no_ident_arb]);
	if ( strchr(buffer,'@') != NULL )
		*strchr(buffer,'@') = 0;
	l = strlen(buffer);
	for ( i = 0 ; i < nb_uniq_species ; ++i )
		if ( strncmp(buffer,uniq_species[i],MIN(l,strlen(uniq_species[i]))) == 0 )
			break;
	if ( i == nb_uniq_species )
		strcpy(uniq_species[nb_uniq_species++],buffer);
	else
		if ( strlen(uniq_species[i]) > l )
			strcpy(uniq_species[i],buffer);
	}
}

/****************************************************************************/
void attribute_clade_uniq_species()
{
int i, j, n;

for ( n = 0 ; n < nb_clade ; ++n )
	nb_uniq_species_per_clade[n] = 0;
for ( i = 0 ; i < nb_uniq_species ; ++i )
	no_clade_uniq_species[i] = -1;

for ( i = 0 ; i < nb_uniq_species ; ++i )
	{
	for ( n = 0 ; n < nb_clade ; ++n )
		{
		for ( j = 0 ; j < nb_species_clade[n] ; ++j )
			{
			if ( strncmp(uniq_species[i],species_clade[n][j],MAX(strlen(uniq_species[i]),strlen(species_clade[n][j]))) == 0 )
				break;
			}
		if ( j != nb_species_clade[n] )
			{
			if ( no_clade_uniq_species[i] == -1 )
				{
				no_clade_uniq_species[i] = n;
				++nb_uniq_species_per_clade[n];
				}
			else
				{
				fprintf(stderr,"%s: %s belongs to several clades\n",file_ARB,uniq_species[i]);
				exit(1);
				}
			}
		}
	}
//for ( i = 0 ; i < nb_uniq_species ; ++i )
//	fprintf(stdout,"%s -> %s\n",uniq_species[i],name_clade[no_clade_uniq_species[i]]);
//for ( i = 0 ; i < nb_clade ; ++i )
//	fprintf(stdout,"%s: %d\n",name_clade[i],nb_uniq_species_per_clade[i]);
}

/****************************************************************************/
void estimate_div_subtree(P_NODE pnode)
{
int i, l;

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	estimate_div_subtree(pnode->gauche);
	estimate_div_subtree(pnode->droite);
	}
else
	{
	strcpy(buffer,ident_arb[pnode->no_ident_arb]);
	if ( strchr(buffer,'@') != NULL )
		*strchr(buffer,'@') = 0;
	l = strlen(buffer);
	for ( i = 0 ; i < nb_uniq_species ; ++i )
		if ( strncmp(buffer,uniq_species[i],l) == 0 )
			break;
	if ( drap_species[i] == -1 )
		{
		drap_species[i] = 0;
		++nb_div_subtree;
		}
	}
}

/****************************************************************************/
void estimate_div_no_subtree(P_NODE pnode, P_NODE psub_tree)
{
int i, l;

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	if ( pnode->gauche != psub_tree )
		estimate_div_no_subtree(pnode->gauche,psub_tree);
	if ( pnode->droite != psub_tree )
		estimate_div_no_subtree(pnode->droite,psub_tree);
	}
else
	{
	strcpy(buffer,ident_arb[pnode->no_ident_arb]);
	if ( strchr(buffer,'@') != NULL )
		*strchr(buffer,'@') = 0;
	l = strlen(buffer);
	for ( i = 0 ; i < nb_uniq_species ; ++i )
		if ( strncmp(buffer,uniq_species[i],l) == 0 )
			break;
	if ( drap_species[i] == -1 )
		{
		drap_species[i] = 0;
		++nb_div_no_subtree;
		}
	}
}

/****************************************************************************/
void find_best_split_div_taxon_int(P_NODE pnode)
{
int i;

for ( i = 0 ; i < nb_uniq_species ; ++i )
	drap_species[i] = -1;
nb_div_subtree = 0;
estimate_div_subtree(pnode);

for ( i = 0 ; i < nb_uniq_species ; ++i )
	drap_species[i] = -1;
nb_div_no_subtree = 0;
estimate_div_no_subtree(proot,pnode);

pnode->nb_taxon_subtree = nb_div_subtree;
pnode->nb_taxon_no_subtree = nb_div_no_subtree;

//fprintf(stdout,"tot=%d Node %d (%f): #subtree=%d #nosubtree=%d\n",nb_div_subtree+nb_div_no_subtree,pnode->no_ident_arb,pnode->dist,nb_div_subtree,nb_div_no_subtree);
if ( max_div <= (nb_div_subtree + nb_div_no_subtree) )
	{
fprintf(stderr,"Node %d (%f): ",pnode->no_ident_arb,pnode->dist);
	if ( max_div == (nb_div_subtree + nb_div_no_subtree) )
		{
fprintf(stderr,"trop court\n");
		if ( pnode->dist < max_internal )
			return;
		}
	pcour = pnode;
	max_div = nb_div_subtree + nb_div_no_subtree;
	max_internal = pnode->dist;
	nb_taxa1 = nb_div_subtree;
	nb_taxa2 = nb_div_no_subtree;
//fprintf(stdout,"MADDIV=%d MAXINTERNAL=%f #taxa1=%d #taxa2=%d (max=%d) (diff_size= %f/%f =%f)\n",max_div,max_internal,nb_taxa1,nb_taxa2,nb_uniq_species,pnode->diff_subtrees,min_diff_subtrees,100.0*pnode->diff_subtrees/min_diff_subtrees);
	}
}

/****************************************************************************/
void find_best_split_div_taxon(P_NODE pnode)
{
find_best_split_div_taxon_int(pnode);

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	find_best_split_div_taxon(pnode->gauche);
	find_best_split_div_taxon(pnode->droite);
	}
}

/****************************************************************************/
void look_for_outoutparalog(P_NODE pnode)
{
int i, n;

for ( i = 0 , n = 0 ; i < nb_clade ; ++i )
	if ( pnode->no_clade_subtree[i] != 0 )
		++n;
//fprintf(stdout,"NODE %d:%d clades\n",pnode->no_ident_arb,n);

//On parcourt l'arbre jusqu'a trouver un clade monophyletique
if ( n == 1 )
	{
//for ( i = 0 ; i < nb_uniq_species ; ++i )	fprintf(stdout,"SP>%s: %d\n",uniq_species[i],proot->no_uniq_species_subtree[i]); fflush(stdout);
// Counts the number of times a species appears in a subtree that corresponds to a monophyletic group, normally one if there are no outoutparalogues (i.e., if the species appears in several subclades, there is outoutparalogy).
	for ( i = 0 ; i < nb_uniq_species ; ++i )
		if ( pnode->no_uniq_species_subtree[i] != 0 )
			++nb_outoutparalog[i];
	}
else
	{
	look_for_outoutparalog(pnode->droite);
	look_for_outoutparalog(pnode->gauche);
	}
}

/****************************************************************************/
void look_for_inparalog(P_NODE pnode)
{
int i, j, n;

for ( i = 0 , n = 0 ; i < nb_uniq_species ; ++i )
	if ( pnode->no_uniq_species_subtree[i] > 0 )
		{
		++n;
		j = i;
		}
//fprintf(stdout,"copies: %d=%d\n",n,pnode->nb_uniq_species_subtree);
if ( n == 1 )
	{
//for ( i = 0 ; i < nb_uniq_species ; ++i )
//	fprintf(stdout,"%2d ",pnode->no_uniq_species_subtree[i]);
//fprintf(stdout,"\n");
	if ( pnode->no_uniq_species_subtree[j] > 1 )
		nb_inparalog[j] += pnode->no_uniq_species_subtree[j];
	}
else
	{
	look_for_inparalog(pnode->droite);
	look_for_inparalog(pnode->gauche);
	}
}

/****************************************************************************/
void compute_statistics_paralogy()
{
int i;

for ( i = 0 , nbtot_outoutparalog = 0 ; i < nb_uniq_species ; ++i )
	if ( nb_outoutparalog[i] > 1 )
		{
		++nbtot_outoutparalog;
		fprintf(stdout,"%s\t%d out-out-paralogs (%s)\n",uniq_species[i],nb_outoutparalog[i],file_ARB);
		}
for ( i = 0 ; i < nb_clade ; ++i)
	no_clade_outout[i] = 0;
for ( i = 0 ; i < nb_uniq_species ; ++i )
	if  ( nb_outoutparalog[i] > 1 )
		{
		++no_clade_outout[no_clade_uniq_species[i]];
/*
		for ( n = 0 ; n < nb_clade ; ++n )
			{
			for ( j = 0 ; j < nb_species_clade[n] ; ++j )
				{
				if ( strncmp(uniq_species[i],species_clade[n][j],MAX(strlen(uniq_species[i]),strlen(species_clade[n][j]))) == 0 )
					break;
				}
			if ( j != nb_species_clade[n] )
				{
				++no_clade_outout[n];
				break;
				}
			}
		if ( n == nb_clade )
			{
			fprintf(stderr,"This species cannot be affected to a clade: %s\n",ident_arb[i]);
			exit(1);
			}
*/
		}
for ( i = 0 , nb_clade_outout = 0 ; i < nb_clade ; ++i)
	if ( no_clade_outout[i] > 0 )
		{
		++nb_clade_outout;
		fprintf(stdout,"Clade %s\t%d out-out-paralogs [%.2f%c] (%s)\n",name_clade[i],no_clade_outout[i],100.0 * (float) no_clade_outout[i] / (float) nb_uniq_species_per_clade[i],'%',file_ARB);
		}

for ( i = 0 ; i < nb_clade ; ++i)
	no_clade_outin[i] = 0;
for ( i = 0 ; i < nb_uniq_species ; ++i )
	if  ( (nb_outinparalog[i] - nb_inparalog[i]) > 1 )
		{
		++no_clade_outin[no_clade_uniq_species[i]];
/*
		for ( n = 0 ; n < nb_clade ; ++n )
			{
			for ( j = 0 ; j < nb_species_clade[n] ; ++j )
				{
				if ( strncmp(uniq_species[i],species_clade[n][j],MAX(strlen(uniq_species[i]),strlen(species_clade[n][j]))) == 0 )
					break;
				}
			if ( j != nb_species_clade[n] )
				{
				++no_clade_outin[n];
				break;
				}
			}
		if ( n == nb_clade )
			{
			fprintf(stderr,"This species cannot be affected to a clade: %s\n",ident_arb[i]);
			exit(1);
			}
*/
		}
for ( i = 0 , nb_clade_outin = 0 ; i < nb_clade ; ++i )
	if ( no_clade_outin[i] > 0 )
		{
		++nb_clade_outin;
		fprintf(stdout,"Clade %s\t%d out-in-paralogs (%s)\n",name_clade[i],no_clade_outin[i],file_ARB);
		}
for ( i = 0 , nbtot_inparalog = 0 ; i < nb_uniq_species ; ++i )
	if ( nb_inparalog[i] > 1 )
		{
		++nbtot_inparalog;
		fprintf(stdout,"%s\t%d in-paralogs (%s)\n",uniq_species[i],nb_inparalog[i],file_ARB);
		}
for ( i = 0 , nbtot_outinparalog = 0 ; i < nb_uniq_species ; ++i )
	{
	if ( nb_inparalog[i] > 1 )
		nb_outinparalog[i] = proot->no_uniq_species_subtree[i] - nb_outoutparalog[i] - nb_inparalog[i] + 1;
	else
		nb_outinparalog[i] = proot->no_uniq_species_subtree[i] - nb_outoutparalog[i] - nb_inparalog[i];
	if ( nb_outinparalog[i] != 0 )
		{
		++nbtot_outinparalog;
		fprintf(stdout,"%s\t%d out-in-paralogs (%s)\n",uniq_species[i],nb_outinparalog[i],file_ARB);
		}
	}

for ( i = 0 , nbtot_outparalog = nb_ident_arb - nb_uniq_species ; i < nb_uniq_species ; ++i )
	if ( nb_inparalog[i] > 1 )
		nbtot_outparalog -= (nb_inparalog[i] - 1);
}

/****************************************************************************/
char is_clade_with_outout_outin_paralogy(int i)
{
if ( (no_clade_outout[i] > 0) || (no_clade_outin[i] > 0) )
	return(0);
else
	return(1);
}

/****************************************************************************/
void create_alignments_paralogy()
{
int i, j;
FILE *out;

system("if ! [ -d ali_paralogy_clade ]; then mkdir ali_paralogy_clade ; fi");

// Creation of the alignment file for clades without paralogy
if ( strrchr(file_ALI,'/') == NULL )
	sprintf(file_SPLIT,"ali_paralogy_clade/%s",file_ALI);
else
	sprintf(file_SPLIT,"ali_paralogy_clade/%s",strrchr(file_ALI,'/') + 1);
*strchr(file_SPLIT,'.') = 0;
sprintf(buffer,"%s-clade-noparalogy.ali",file_SPLIT);
out = ffopen(buffer,"w");
fprintf(out,"#%s\n#%s\n",molecule,domaine);
for ( i = 0 ; i < nb_ident_ali ; ++i )
	if ( is_clade_with_outout_outin_paralogy(get_no_clade_ident(ident_ali + i * LENIDENT)) == 1 )
		fprintf(out,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
fclose(out);
for ( j = 0 ; j < nb_clade ; ++j )
	if ( is_clade_with_outout_outin_paralogy(j) == 0 )
		{
		sprintf(buffer,"%s-clade-%s.ali",file_SPLIT,name_clade[j]);
		out = ffopen(buffer,"w");
		fprintf(out,"#%s\n#%s\n",molecule,domaine);
		for ( i = 0 ; i < nb_ident_ali ; ++i )
			if ( get_no_clade_ident(ident_ali + i * LENIDENT) == j )
				fprintf(out,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
		fclose(out);
		}
}


/******************************************************************************/
/* Function affecting the taxonomic distribution of each node                 */
/******************************************************************************/
void compute_taxonomy_complement_subtree(P_NODE pnode)
{
int i;

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	for ( i = 0 ; i < nb_uniq_species ; ++i )
		{
		(pnode->gauche)->no_uniq_species_complement_subtree[i] = pnode->no_uniq_species_complement_subtree[i] + (pnode->droite)->no_uniq_species_subtree[i];
		(pnode->droite)->no_uniq_species_complement_subtree[i] = pnode->no_uniq_species_complement_subtree[i] + (pnode->gauche)->no_uniq_species_subtree[i];
		}
	compute_taxonomy_complement_subtree(pnode->gauche);
	compute_taxonomy_complement_subtree(pnode->droite);
	}
}

/****************************************************************************/
void look_for_clan_int(P_NODE pnode)
{
int i, n;
int nb_uniq_species_per_node[MAXSPECIES];

//Tests one side of the bipartition
for ( i = 0 ; i < nb_clade ; ++i )
	nb_uniq_species_per_node[i] = 0;
for ( i = 0 ; i < nb_uniq_species ; ++i )
	nb_uniq_species_per_node[no_clade_uniq_species[i]] += pnode->no_uniq_species_subtree[i];

for ( n = 0 ; n < nb_clade ; ++n )
	if ( (nb_uniq_species_per_clan[n] != 0) && (drap_presence_clan[n] == 0) )
		{
		if ( nb_uniq_species_per_clan[n] == nb_uniq_species_per_node[n] )
			{
			for ( i = 0 ; i < nb_clade ; ++i )
				if ( (i != n) && (nb_uniq_species_per_node[i] != 0) )
					break;
			if ( i == nb_clade )
				{
				drap_presence_clan[n] = 1;
/*
fprintf(stdout,"\nNoeud=%3d Dist=%.2f Size=%.2f\n\t",pnode->no_ident_arb,pnode->dist,pnode->size_subtree); fflush(stdout);
for ( i = 0 ; i < nb_uniq_species ; ++i )
	fprintf(stdout,"%d-",pnode->no_uniq_species_subtree[i]);
fprintf(stdout,"\n\t");
for ( i = 0 ; i < nb_uniq_species ; ++i )
	fprintf(stdout,"%d-",pnode->no_uniq_species_complement_subtree[i]);
fprintf(stdout,"\n");
for ( i = 0 ; i < nb_clade ; ++i )
	fprintf(stdout,"%s: %d unique species\n",name_clade[i],nb_uniq_species_per_node[i]);
*/
				if ( (nb_uniq_species_per_clan[n] <= 1) || ((nb_ident_arb - nb_uniq_species_per_clan[n]) <= 1) )
					fprintf(stdout,"%s: CLAN %s trivially present\n",file_ARB,name_clade[n]);
				else
					fprintf(stdout,"%s: CLAN %s present\n",file_ARB,name_clade[n]);
				}
			}
		}

//Tests the other side of the bipartition
for ( i = 0 ; i < nb_clade ; ++i )
	nb_uniq_species_per_node[i] = 0;
for ( i = 0 ; i < nb_uniq_species ; ++i )
	nb_uniq_species_per_node[no_clade_uniq_species[i]] += pnode->no_uniq_species_complement_subtree[i];

for ( n = 0 ; n < nb_clade ; ++n )
	if ( (nb_uniq_species_per_clan[n] != 0) && (drap_presence_clan[n] == 0) )
		{
		if ( nb_uniq_species_per_clan[n] == nb_uniq_species_per_node[n] )
			{
			for ( i = 0 ; i < nb_clade ; ++i )
				if ( (i != n) && (nb_uniq_species_per_node[i] != 0) )
					break;
			if ( i == nb_clade )
				{
				drap_presence_clan[n] = 1;
/*
fprintf(stdout,"\nNoeud=%3d Dist=%.2f Size=%.2f\n\t",pnode->no_ident_arb,pnode->dist,pnode->size_subtree); fflush(stdout);
for ( i = 0 ; i < nb_uniq_species ; ++i )
	fprintf(stdout,"%d-",pnode->no_uniq_species_subtree[i]);
fprintf(stdout,"\n\t");
for ( i = 0 ; i < nb_uniq_species ; ++i )
	fprintf(stdout,"%d-",pnode->no_uniq_species_complement_subtree[i]);
fprintf(stdout,"\n");
for ( i = 0 ; i < nb_clade ; ++i )
	fprintf(stdout,"%s: %d unique species\n",name_clade[i],nb_uniq_species_per_node[i]);
*/
				if ( (nb_uniq_species_per_clan[n] <= 1) || ((nb_ident_arb - nb_uniq_species_per_clan[n]) <= 1) )
					fprintf(stdout,"%s: CLAN %s trivially present\n",file_ARB,name_clade[n]);
				else
					fprintf(stdout,"%s: CLAN %s present\n",file_ARB,name_clade[n]);
				}
			}
		}

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	look_for_clan_int(pnode->gauche);
	look_for_clan_int(pnode->droite);
	}
}

/****************************************************************************/
void look_for_clan()
{
int i, n;

for ( i = 0 ; i < nb_uniq_species ; ++i )
	{
	(proot->gauche)->no_uniq_species_complement_subtree[i] = (proot->droite)->no_uniq_species_subtree[i];
	(proot->droite)->no_uniq_species_complement_subtree[i] = (proot->gauche)->no_uniq_species_subtree[i];
	}
compute_taxonomy_complement_subtree(proot->gauche);
compute_taxonomy_complement_subtree(proot->droite);

for ( i = 0 ; i < nb_clade ; ++i )
	{
	drap_presence_clan[i] = 0;
	nb_uniq_species_per_clan[i] = 0;
	}
for ( i = 0 ; i < nb_uniq_species ; ++i )
	nb_uniq_species_per_clan[no_clade_uniq_species[i]] += proot->no_uniq_species_subtree[i];
//for ( i = 0 ; i < nb_clade ; ++i )
//	fprintf(stdout,"%s: %d unique species\n",name_clade[i],nb_uniq_species_per_clan[i]);
look_for_clan_int(proot);
for ( i = 0 , n = 0 ; i < nb_clade ; ++i )
	if ( nb_uniq_species_per_clan[i] > 1 )
		if ( drap_presence_clan[i] == 0 )
			{
			++n;
			fprintf(stdout,"%s: CLAN %s absent\n",file_ARB,name_clade[i]);
			}

}

/****************************************************************************/
void compte_nb_species(P_NODE pnode, int *nb_seq)
{
if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	compte_nb_species(pnode->gauche,nb_seq);
	compte_nb_species(pnode->droite,nb_seq);
	}
else
	++*nb_seq;
}

/****************************************************************************/
void display_species_long_internal_branch(P_NODE pnode, P_NODE pref, char drap_complementary)
{
if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	if ( (drap_complementary == 1) && (pnode == pref) )
		return;
	display_species_long_internal_branch(pnode->gauche,pref,drap_complementary);
	display_species_long_internal_branch(pnode->droite,pref,drap_complementary);
	}
else
	fprintf(stdout,"remove_seq \"%s\" %s\t#internal BRANCH too long %.2f\n",ident_arb[pnode->no_ident_arb],file_ALI,pref->dist);
}

/****************************************************************************/
void look_for_long_branch(P_NODE pnode)
{
int n;

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	look_for_long_branch(pnode->gauche);
	look_for_long_branch(pnode->droite);
	if ( pnode->dist > min_long_BL )
		{
//sprintf(file_ARB,"A.arb");
//save_tree_arb();
		n = 0;
		compte_nb_species(pnode,&n);
		fprintf(stdout,"%s: too long internal BRANCH (%d - %d species) %.2f\n",file_ALI,pnode->no_ident_arb,MIN(n,nb_ident_arb - n),pnode->dist);
		if ( n < (nb_ident_arb / 2) )
			display_species_long_internal_branch(pnode,pnode,0);
		else
			display_species_long_internal_branch(proot,pnode,1);
		}
	}
else
	if ( pnode->dist > min_long_BL )
		{
		fprintf(stdout,"remove_seq \"%s\" %s\t#too long terminal branch %.2f\n",ident_arb[pnode->no_ident_arb],file_ALI,pnode->dist);
		}
}

/****************************************************************************/
void usage(void)
{
fprintf(stderr,"USAGE: orthology-check-ali");
fprintf(stderr,"\n\tarb=<name of file .arb>");
fprintf(stderr,"\n\tali=<name of file .ali>");
fprintf(stderr,"\n\tclades=<file containing the species list for each clade>");
fprintf(stderr,"\n\tmaxdist_sistergroup=<maximum distance to consider sister-group relationship as suspicious (default=0))>");
fprintf(stderr,"\n\tmindist_suspicious=<minimum distance to consider an internal branch length as sufficient to support suspicious taxonomy (default=0)>");
fprintf(stderr,"\n\tmin_dist_subclade=<maximum distance to consider 2 subgroups of a clade as putative paralogs (default=0)>");
fprintf(stderr,"\n\tmin_nb_spec_split=<minimum number of species in the smallest subclades to allow splitting (default=1)>");
fprintf(stderr,"\n\tmin_overlap=<number of amino acids overlapping between two subclades to validate a putative paralogy (default=50)>\n\t");
fprintf(stderr,"\n\tmin_long_BL=<minimum relative branch length to consider a branch as too long (default=50)>\n\t");
fprintf(stderr,"\n\tparalogy_clade=<false or true> (to create separate alignments for clades without out- and in-paralogy and for clades with out- or in-paralogy)\n\n");
exit(1);
}

/****************************************************************************/
int main(int argc, char **argv)
{
int i, j, k,  nb0, nb1, nb2, n1, n2, nbtot_subclade, max_subclade;
FILE *out1, *out2, *out_db;
float d, max_d, best1, best2;
char drap_paralogy[MAXCLADE], *ps, drap_splitting, drap_paralogy_clade; //drap_splitting = 0 (no split) 1 (split with the longest branch when n>2 subclades exist) 2 (when the clades are at most splitted into 2 subclades and all the clades share the same longest branch)))

maxdist_sistergroup = 0.0;
mindist_suspicious = 0.0;
min_dist_subclade = 0.0;
min_overlap = 50;
min_nb_spec_split = 1;
drap_paralogy_clade = 0;	//indicating that no alignments related to paralogy will be created
min_long_BL = 50.0;
for( i = 1 ; i < argc ; ++i )
	{
	if ( strncmp(argv[i],"ali=",4) == 0 )
		{
		strncpy(file_ALI,argv[i] + 4,LENFILENAME - 1);
		if ( strchr(file_ALI,'.') != NULL )
			*strchr(file_ALI,'.') = 0;
		strcat(file_ALI,".ali");
		}
	else
		if ( strncmp(argv[i],"arb=",4) == 0 )
			{
			strncpy(file_ARB,argv[i] + 4,LENFILENAME - 1);
			if ( strchr(file_ARB,'.') != NULL )
				*strchr(file_ARB,'.') = 0;
			strcat(file_ARB,".arb");
			}
		else
			if ( strncmp(argv[i],"clades=",7) == 0 )
				strcpy(file_clade,argv[i] + 7);
			else
				if ( strncmp(argv[i],"maxdist_sistergroup=",20) == 0 )
					sscanf(argv[i],"maxdist_sistergroup=%f",&maxdist_sistergroup);
				else
					if ( strncmp(argv[i],"mindist_suspicious=",19) == 0 )
						sscanf(argv[i],"mindist_suspicious=%f",&mindist_suspicious);
					else
						if ( strncmp(argv[i],"min_dist_subclade=",18) == 0 )
							sscanf(argv[i],"min_dist_subclade=%f",&min_dist_subclade);
						else
							if ( strncmp(argv[i],"min_overlap=",12) == 0 )
								sscanf(argv[i],"min_overlap=%d",&min_overlap);
							else
								if ( strncmp(argv[i],"min_nb_spec_split=",17) == 0 )
									sscanf(argv[i],"min_nb_spec_split=%d",&min_nb_spec_split);
								else
									if ( strncmp(argv[i],"paralogy_clade=false",20) == 0 )
										drap_paralogy_clade = 0;
									else
										if ( strncmp(argv[i],"paralogy_clade=true",19) == 0 )
											drap_paralogy_clade = 1;
										else
											if ( strncmp(argv[i],"min_long_BL=",12) == 0 )
												sscanf(argv[i],"min_long_BL=%f",&min_long_BL);
											else
												{
												fprintf(stderr,"This argument is not correctly formatted: %s\n\n",argv[i]);
												usage();
												}
	}

//fprintf(stdout,"OK-4\n"); fflush(stdout);
read_file_CLADE();
//fprintf(stdout,"OK-3\n"); fflush(stdout);
read_file_ALI();
//fprintf(stdout,"OK-2\n"); fflush(stdout);
if ( read_File_ARB() == -1 )
	exit(-1);
//fprintf(stdout,"OK-1\n"); fflush(stdout);
allocate_no_clade();
//fprintf(stdout,"OK0\n"); fflush(stdout);

if ( (proot = compute_tree()) == NULL )
	exit(-1);
// Roots on species 0, with one of its two branches at 0, to ensure the midlength search
//fprintf(stdout,"OK1\n"); fflush(stdout);
root_left_species0();
compute_size_subtree(proot);
//fprintf(stdout,"OK2\n"); fflush(stdout);
// Searches for midpoint rooting and stores the size differences between the two subtrees
min_diff_subtrees = proot->size_subtree;
pcour = proot;
lookfor_midlength(proot);
percent_BL_root = 0.5;
root_left();
compute_size_subtree(proot);
normalise100_tree(proot,proot->size_subtree);
//fprintf(stdout,"OK3\n"); fflush(stdout);

//verif_tree(proot);
//sprintf(file_ARB,"A.arb");
//save_tree_arb();

// to compute the number of out-out-paralogs and of inparalogs
nb_uniq_species = 0;
generate_liste_uniq_species(proot);
attribute_clade_uniq_species();
compute_taxonomy_subtree(proot);
//fprintf(stdout,"OK6\n"); fflush(stdout);
//for ( i = 0 ; i < nb_uniq_species ; ++i )	fprintf(stdout,"SP>%s: %d\n",uniq_species[i],proot->no_uniq_species_subtree[i]); fflush(stdout);
//for ( i = 0 ; i < nb_uniq_species ; ++i )
//	drap_species[i] = -1;
//nb_div_subtree = 0;
//estimate_div_subtree(pnode);
//verif_tree(proot);
look_for_clan();
for ( i = 0 ; i < nb_uniq_species ; ++i )
	{
	nb_outoutparalog[i] = 0;
	nb_inparalog[i] = 0;
	nb_outinparalog[i] = proot->no_uniq_species_subtree[i];
	}
look_for_outoutparalog(proot);
look_for_inparalog(proot);
compute_statistics_paralogy();
if ( drap_paralogy_clade == 1 )
	create_alignments_paralogy();

lookfor_suspicious_sistergroup(proot);
nb_monophyletic_clade = 0;
lookfor_clade_suspicious_taxonomy(proot);
if ( nb_monophyletic_clade != nb_nonempty_clade)
	{
	max_subclade = -1;
	drap_splitting = 0;
// Looking for cases of non-monophyly that fulfill the criteria of putative paralogy
	for ( i = 0 , nbtot_subclade = 0; i < nb_clade ; ++i )
		{
		drap_paralogy[i] = -1;
		if ( (drap_monophyly_clade[i] == -1) && (no_clade_tree[i] > 1) )
			{
			nb_monophyletic_subclade = 0;
			for ( j = 0 ; j < MAXSUBCLADE ; ++j )
				p_subclade[j] = NULL;
			lookfor_monophyletic_subclade(proot,i);
			fprintf(stdout,"%s: %d subclades\n",name_clade[i],nb_monophyletic_subclade);
			nbtot_subclade += nb_monophyletic_subclade;
			for ( j = 0 , max_d = 0 ; j < nb_monophyletic_subclade ; ++j )
				for ( k = 0 ; k < j ; ++k )
					{
					compute_dist_of_two_nodes(p_subclade[j],p_subclade[k],&d,i);
					if ( lookfor_overlap_subclades(p_subclade[j],p_subclade[k]) == 0 )
						max_d = MAX (max_d,d);
					}
//fprintf(stdout,"%s: LCA12=%d\n",name_clade[i],p_LCA12[i]->no_ident_arb);
			if ( max_d >= min_dist_subclade )
				{
				drap_paralogy[i] = nb_monophyletic_subclade;
				max_subclade = MAX(max_subclade,nb_monophyletic_subclade);
				fprintf(stdout,"%s: putative paralogy for clade %s (max_dist=%.2f)\n",file_ARB,name_clade[i],max_d);
				}
			}
		}
//fprintf(stdout,"max=%d\n",max_subclade);
	fprintf(stdout,"%s: %d monophyletic clades over %d (a total of %d subclades) [%d species with inparalogs, %d out-in-paralogs (involving %d clades) and %d out-out-paralogs (involving %d clades)]{%d out-paralogs (%.0f%c)}\n",file_ARB,nb_monophyletic_clade,nb_nonempty_clade,nbtot_subclade,nbtot_inparalog,nbtot_outinparalog,nb_clade_outin,nbtot_outoutparalog,nb_clade_outout,nbtot_outparalog,100.0 * (float) nbtot_outparalog / (float) nb_uniq_species,'%');

//for ( i = 0 ; i < nb_clade ; ++i )
//	if ( drap_paralogy[i] != -1 )
//		fprintf(stdout,"%slongest_branch=%d (%.2f)\n",name_clade[i],p_longest_branch[i]->no_ident_arb,p_longest_branch[i]->dist);

	if ( max_subclade == 2 )
		{
// Checking that the longest branch is the same for all the clades
		drap_splitting = 2;
		for ( i = 0 , pcour = NULL ; i < nb_clade ; ++i )
			if ( drap_paralogy[i] != -1 )
				{
				if ( pcour == NULL )
					pcour = p_longest_branch[i];
				else
					{
					if ( pcour != p_longest_branch[i] )
						{
						fprintf(stdout,"File %s has a MAJOR PARALOGY PROBLEM since the long branch is not the same for the different clades\n",file_ARB);
						drap_splitting = 0;
						break;
						}
					}
				}
		root_left();
		}
	else
		if ( max_subclade > 2 )
			{
// Looking for the longest branch for all the pairs of subclades
			drap_splitting = 1;
			pcour = NULL;
			for ( i = 0 ; i < nb_clade ; ++i )
				{
				if ( (drap_monophyly_clade[i] == -1) && (no_clade_tree[i] > 1) )
					{
					nb_monophyletic_subclade = 0;
					for ( j = 0 ; j < MAXSUBCLADE ; ++j )
						p_subclade[j] = NULL;
					lookfor_monophyletic_subclade(proot,i);
					for ( j = 0 ; j < nb_monophyletic_subclade ; ++j )
						for ( k = 0 ; k < j ; ++k )
							{
							compute_dist_of_two_nodes(p_subclade[j],p_subclade[k],&d,i);
//fprintf(stdout,"%s: longest=%d (%.2f)\n",name_clade[i],p_longest_branch[i]->no_ident_arb,p_longest_branch[i]->dist);
							if ( lookfor_overlap_subclades(p_subclade[j],p_subclade[k]) == 0 )
								{
								if ( pcour == NULL )
									pcour = p_longest_branch[i];
								else
									{
									if ( p_longest_branch[i]->dist > pcour->dist )
										pcour = p_longest_branch[i];
									}
								}
							}
					}
				}
//fprintf(stdout,"%s: longest overall=%d (%.2f)\n",name_clade[i],pcour->no_ident_arb,pcour->dist);
			root_left();
			}

// Splitting of the alignment if possible
	strcpy(buffer,file_ARB);
	if ( drap_splitting == 1 )
		{
		system("if ! [ -d SPLITTEDl ]; then mkdir SPLITTEDl ; fi");
		sprintf(file_ARB,"SPLITTEDl/%s",buffer);
		save_tree_arb();
		}
	else
		if ( drap_splitting == 2 )
			{
			system("if ! [ -d SPLITTED2 ]; then mkdir SPLITTED2 ; fi");
			sprintf(file_ARB,"SPLITTED2/%s",buffer);
			save_tree_arb();
			}

// Cases where one can split
// One computes the number of species in each subtree
	erase_taxonomy_subtree(proot);
	compute_taxonomy_subtree(proot);
	for ( i = 0 , n1 = 0 , n2 = 0 ; i < nb_uniq_species ; ++i )
		{
		if ( (proot->gauche)->no_uniq_species_subtree[i] >= 1 )
			++n1;
		if ( (proot->droite)->no_uniq_species_subtree[i] >= 1 )
			++n2;
//fprintf(stdout,"%d ",(proot->droite)->no_uniq_species_subtree[i]);
		}
//fprintf(stdout,"\n");
	if ( (drap_splitting > 0) && (((proot->gauche)->dist + (proot->droite)->dist) >= min_dist_subclade) && (n1 >= min_nb_spec_split) && (n2 >= min_nb_spec_split) )
		{
//fprintf(stdout,"(proot->gauche)->no_uniq_species_subtree=%d et (proot->droite)->no_uniq_species_subtree=%d\n",n1,n2);
//fprintf(stdout,"(proot->gauche)->nb_taxon_subtree=%d et (proot->droite)->nb_taxon_subtree=%d\n",(proot->gauche)->nb_taxon_subtree,(proot->droite)->nb_taxon_subtree);
// to identify in which file .ali the sequences should be placed
		nb0 = nb1 = nb2 = 0;
		nb_spec_subtree = 0;
		lookfor_species_subtree(proot->gauche);
		for ( i = 0 ; i < nb_ident_ali ; ++i )
			{
			drap_species[i] = 0;
			for ( j = 0 ; j < nb_spec_subtree ; ++j )
				if ( strncmp(ident_ali + i * LENIDENT,ident_arb[no_spec_subtree[j]],strlen(ident_arb[no_spec_subtree[j]])) == 0 )
					break;
			if ( j != nb_spec_subtree )
				{
				drap_species[i] = 1;
				++nb1;
				}
			}
		nb_spec_subtree = 0;
		lookfor_species_subtree(proot->droite);
		for ( i = 0 ; i < nb_ident_ali ; ++i )
			{
			for ( j = 0 ; j < nb_spec_subtree ; ++j )
				if ( strncmp(ident_ali + i * LENIDENT,ident_arb[no_spec_subtree[j]],strlen(ident_arb[no_spec_subtree[j]])) == 0 )
					break;
			if ( j != nb_spec_subtree )
				{
				drap_species[i] = 2;
				++nb2;
				}
			}

		if ( drap_splitting == 1 )
			{
			if ( strrchr(file_ALI,'/') == NULL )
				sprintf(file_SPLIT,"SPLITTEDl/%s",file_ALI);
			else
				sprintf(file_SPLIT,"SPLITTEDl/%s",strrchr(file_ALI,'/') + 1);
			}
		else
			{
			if ( strrchr(file_ALI,'/') == NULL )
				sprintf(file_SPLIT,"SPLITTED2/%s",file_ALI);
			else
				sprintf(file_SPLIT,"SPLITTED2/%s",strrchr(file_ALI,'/') + 1);
			}
		*strchr(file_SPLIT,'.') = 0;
		sprintf(buffer,"%s-1.ali",file_SPLIT);
		out1 = ffopen(buffer,"w");
		sprintf(buffer,"%s-2.ali",file_SPLIT);
		out2 = ffopen(buffer,"w");
		fprintf(out1,"#%s\n#%s\n",molecule,domaine);
		fprintf(out2,"#%s\n#%s\n",molecule,domaine);
		out_db = ffopen("db_tmp","w");
		for ( i = 0 ; i < nb_ident_ali ; ++i )
			if ( drap_species[i] == 1 )
				{
				fprintf(out1,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
				fprintf (out_db, ">db1|%d\n", i);
				for ( ps = seq + i * lenseq_tot ; *ps != 0 ; ++ps )
					if ( strchr("*-$ X?",*ps) == NULL )
						fprintf (out_db, "%c",*ps);
				fprintf (out_db, "\n");
				}
			else
				if ( drap_species[i] == 2 )
					{
					fprintf(out2,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
					fprintf (out_db, ">db2|%d\n", i);
					for ( ps = seq + i * lenseq_tot ; *ps != 0 ; ++ps )
						if ( strchr("*-$ X?",*ps) == NULL )
							fprintf (out_db, "%c",*ps);
					fprintf (out_db, "\n");
					}
		fclose(out1);
		fclose(out2);
		fclose(out_db);

//putting the SEQUENCES OF THE .ALI that are not present in the .arbre
		for ( i = 0 ; i < nb_ident_ali ; ++i )
			if ( drap_species[i] == 0 )
				{
				if ( nb0++ == 0 )
					system("formatdb -i db_tmp -p T -o F");
				out1 = ffopen("seed.tmp","w");
				fprintf(out1, ">seed%d\n", i);
				for ( ps = seq + i * lenseq_tot ; *ps != 0 ; ++ps )
					if ( strchr("*-$ X?",*ps) == NULL )
						fprintf (out1, "%c",*ps);
				fprintf (out1, "\n");
				fclose(out1);
				system("blastp -db db_tmp -query seed.tmp -out blastp.tmp -seg no");
				out1 = ffopen("blastp.tmp","r");
				best1 = best2 = -1.0;
				while ( fgets(buffer,LENBUFFER - 1, out1) != NULL )
					{
					if ( best1 == -1.0 )
						if ( strncmp(buffer,"  db1",5) == 0 )
							{
							sscanf(buffer + 50,"%f",&best1);
							if ( best2 != -1.0 )
								break;
							}
					if ( best2 == -1.0 )
						if ( strncmp(buffer,"  db2",5) == 0 )
							{
							sscanf(buffer + 50,"%f",&best2);
							if ( best1 != -1.0 )
								break;
							}
					}
				fclose(out1);
				if ( best1 > best2 )
					{
					sprintf(buffer,"%s-1.ali",file_SPLIT);
					out1 = ffopen(buffer,"a");
					fprintf(out1,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
					fclose(out1);
					++nb1;
//fprintf(stdout,"sequence absent from .arb and put in -1.ali:%s\n",ident_ali + i * LENIDENT);
					}
				else
					{
					sprintf(buffer,"%s-2.ali",file_SPLIT);
					out2 = ffopen(buffer,"a");
					fprintf(out2,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
					fclose(out2);
					++nb2;
//fprintf(stdout,"sequence absent from .arb and put in -2.ali:%s\n",ident_ali + i * LENIDENT);
					}
				}
		if ( nb1 < nb2 )
			{
			sprintf(buffer,"mv %s-1.ali detect-problems-arb.TMP",file_SPLIT);
			system(buffer);
			sprintf(buffer,"mv %s-2.ali %s-1.ali",file_SPLIT,file_SPLIT);
			system(buffer);
			sprintf(buffer,"mv detect-problems-arb.TMP %s-2.ali",file_SPLIT);
			system(buffer);
			i = nb1;
			nb1 = nb2;
			nb2 = i;
			}
		fprintf(stdout,"File %s has been SPLITTED (%d-%d species)\n",file_ALI,nb1,nb2);
		}
	else
//		fprintf(stdout,"File %s: the internal branch is too short for being splitted (%.2f)\n",file_ALI,(proot->gauche)->dist + (proot->droite)->dist);
		fprintf(stdout,"File %s: the internal branch is too short for being splitted (%.2f < %.2f) or a MAJOR paralogy issue\n",file_ALI,(proot->gauche)->dist + (proot->droite)->dist,min_dist_subclade);
	}
else
	fprintf(stdout,"%s: ALL clades are monophyletic (%d) [%d species with inparalogs, %d out-in-paralogs (involving %d clades) and %d out-out-paralogs (involving %d clades)] {%d out-paralogs (%.0f%c)}\n",file_ARB,nb_monophyletic_clade,nbtot_inparalog,nbtot_outinparalog,nb_clade_outin,nbtot_outoutparalog,nb_clade_outout,nbtot_outparalog,100.0 * (float) nbtot_outparalog / (float) nb_uniq_species,'%');

//	fprintf(stdout,"%s: %d monophyletic clades over %d (a total of %d subclades) [%d species with inparalogs, %d out-in-paralogs and %d out-out-paralogs (involving %d clades)]\n",file_ARB,nb_monophyletic_clade,nb_nonempty_clade,nbtot_subclade,nbtot_inparalog,nbtot_outinparalog,nbtot_outoutparalog,nb_clade_outout);


if ( (proot->gauche)->gauche == NULL )
	{
	(proot->gauche)->dist += (proot->droite)->dist;
	(proot->droite)->dist = 0.0;
	}
else
	{
	(proot->droite)->dist += (proot->gauche)->dist;
	(proot->gauche)->dist = 0.0;
	}
look_for_long_branch(proot);

exit(0);
/*
nb_internal_branch = 0;
lookfor_internal_distance(proot,&nb_internal_branch);
qsort(internal_length,(size_t) nb_internal_branch,sizeof(float),(const void *) comparer_float);
for ( i = 0 ; i < nb_internal_branch ; ++i )
	fprintf(stderr,"%.3f ",internal_length[i]);
fprintf(stderr,"--> %.3f\n",internal_length[(int) (LIMINTERNAL * (float) nb_internal_branch)]);


//Finds the node that maximizes taxonomic diversity on each side of its branch
nb_uniq_species = 0;
max_internal = 0.0;
generate_liste_species(proot);
//for ( i = 0 ; i < nb_uniq_species ; ++i )	fprintf(stdout,"%s\n",uniq_species[i]);
pcour = proot;
max_div = 0;
nb_taxa1 = nb_taxa2 = 0;
find_best_split_div_taxon(proot->gauche);
find_best_split_div_taxon(proot->droite);

if ( ((nb_taxa1 == nb_uniq_species) && (nb_taxa2 > 1)) || ((nb_taxa2 == nb_uniq_species) && (nb_taxa1 > 1)) )
	drap_split = 3;			// cases where at least one subgroup has all taxa
else
	if ( (nb_taxa1 > 10) && (nb_taxa2 > 10) )
		drap_split = 2;			// cases where both subgroups have enough taxa
	else
		drap_split = 0;			// cases where we don't split

//if ( drap_split != 0 )
//	{
	percent_BL_root = 0.5;
	root_left();
//fprintf(stdout,"APRES %d %f\n",pcour->no_ident_arb,pcour->diff_subtrees);
//	}

if ( drap_split != 0 )
	if ( max_internal < internal_length[(int) (LIMINTERNAL * (float) nb_internal_branch)] )
		{
		fprintf(stdout,"Splittable, but too short internal branch\t");
		drap_split = 1;
		}

if ( drap_split <= 1 )
	{
	if ( drap_split == 0 )
		fprintf(stdout,"File %s has not been splitted\n",file_ALI);
	else
		fprintf(stdout,"The splittable file %s has not been splitted\n",file_ALI);
	if ( drap_split == 0 )
		sprintf(buffer,"cp %s unsplitted",file_ALI);
	else
		sprintf(buffer,"cp %s splittable",file_ALI);
	system(buffer);
	strcpy(buffer,file_ARB);
	if ( drap_split == 0 )
		sprintf(file_ARB,"unsplitted/%s",buffer);
	else
		sprintf(file_ARB,"splittable/%s",buffer);
	save_tree_arb();
	}
else
	{
// These two lines allow for good values of nb1 and nb2
	lookfor_species_subtree(proot->gauche);
	lookfor_dup_species(proot->droite);
	nb1 = nb2 = 0;
	fprintf(stdout,"File %s has been splitted",file_ALI);
	if ( drap_split == 2 )
		sprintf(file_SPLIT,"splitted/%s",file_ALI);
	else
		sprintf(file_SPLIT,"complete/%s",file_ALI);
	*strchr(file_SPLIT,'.') = 0;
	sprintf(buffer,"%s-1.ali",file_SPLIT);
	out1 = ffopen(buffer,"w");
	sprintf(buffer,"%s-2.ali",file_SPLIT);
	out2 = ffopen(buffer,"w");
	fprintf(out1,"#%s\n#%s\n",molecule,domaine);
	fprintf(out2,"#%s\n#%s\n",molecule,domaine);
	for ( i = 0 ; i < nb_ident_ali ; ++i )
		{
		for ( j = 0 ; j < nb_spec_subtree ; ++j )
			if ( strcmp(ident_ali + i * LENIDENT,ident_arb[no_spec_subtree[j]]) == 0 )
				break;
		if ( j == nb_spec_subtree )
			{
			fprintf(out2,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
			++nb2;
			}
		else
			{
			fprintf(out1,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
//fprintf(stdout,"%s\n",ident_arb[no_spec_subtree[j]]);
			++nb1;
			}
		}
	fclose(out1);
	fclose(out2);
	strcpy(buffer,file_ARB);
	if ( drap_split == 2 )
		sprintf(file_ARB,"splitted/%s",buffer);
	else
		sprintf(file_ARB,"complete/%s",buffer);
	save_tree_arb();
	fprintf(stdout," (%d and %d species)\t",nb1,nb2);
	if ( drap_split == 2 )
		fprintf(stdout,"incomplete\t");
	else
		fprintf(stdout,"complete\t");
fprintf(stdout,"Node %d: MADDIV= %d #taxa1= %d #taxa2= %d (max= %d) (diff_size= %f/%f %.0f) (internal= %f /%f)\n",pcour->no_ident_arb,max_div,pcour->nb_taxon_subtree,pcour->nb_taxon_no_subtree,nb_uniq_species,pcour->diff_subtrees,min_diff_subtrees,100.0*pcour->diff_subtrees/min_diff_subtrees,pcour->dist,internal_length[(int) (LIMINTERNAL * (float) nb_internal_branch)]);
	}
*/
exit(0);
}
