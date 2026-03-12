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
#define MAXSPECIE		10000
#define LENARBRE		2000000
#define LENBUFFER		LENARBRE
#define LENSEQ			10000000
#define MAXCOMMENT		50    /* maximum number of comment lines */
#define LENCOMMENT		100
#define LENFILENAME		500
#define MAXCARSEQ		500000000
#define LENMOLDOM		100
#define MAXCRITERIA		100
//#define percent_longer_internal_branch		0.10	// percentage of internal branches that can be longer than the re-rooting one

typedef struct node{        /* structure defining a node  */
    float dist;
    float size_subtree;
	float diff_subtrees;
    int no_ident;
	int nb_taxon_subtree;
	int nb_taxon_no_subtree;
	int div_taxon;
    struct node *pere;
    struct node *gauche;
    struct node *droite;
} NODE, *P_NODE;

char file_ARB[LENFILENAME], file_ALI[LENFILENAME], file_SPLIT[LENFILENAME];
int nb_comment;
char comment[MAXCOMMENT][LENCOMMENT];
char arbre[LENARBRE];
char buffer[LENBUFFER];
int enracine;      /* is 1 if the tree is rooted */
int nb_node;
int nb_neg;
char *pc_recurs;
int nb_multif;
int nb_spec;
float percent, min_diff_subtrees;
char ident_arb[MAXSPECIE][LENIDENT];    /* identifier table  */
char species[MAXSPECIE][LENIDENT];    /* identifier table  */
char drap_species[MAXSPECIE];
int nb_species_taxon, nb_div_subtree, nb_div_no_subtree;
int max_div, nb_taxa1, nb_taxa2;
int nb_spec_subtree;
int no_spec_subtree[MAXSPECIE];
int drap_dup, drap_split;
P_NODE proot;            /* pointer to the root structure */
P_NODE pcour;            /* pointer to the current node on the screen  */
int nb_internal_branch;
float internal_BL[MAXSPECIE], max_internal_BL;		// max_internal_BL contains the length of the longest branch that maximizes taxonomic diversity
int nb_criteria, min_nb_species_group1[MAXCRITERIA], min_nb_species_group2[MAXCRITERIA], min_shared_taxa[MAXCRITERIA];
float percent_longer_internal_branch[MAXCRITERIA];
FILE *out_log;

char molecule[LENMOLDOM], domaine[LENMOLDOM];
int type_seq_ali;
int nb_spec_ali;
int nb_nuc_tot;
long lenseq_tot;           /* Maximum length that sequences can have */
                           /* due to MAXCARSEQ and nb_spec_ali */
char seq[MAXCARSEQ];
char *ident_ali;
int mask[MAXSPECIE];

/* Definition of TREEPLOT functions */
void cree_dichotomie(void);
void parse(int level);
char traite_dist(char **pb, char **pa);
char traite_ident(char **pb, char **pa);
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
void save_tree_arb(char *name_DIR);
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
/* Opens a file and verifies that it opened successfully                  */
/******************************************************************************/
FILE *ffopen(char *File, char *Mode)
{
FILE *in;

if ( (in = fopen(File,Mode)) == NULL )
	{
	fprintf(stderr,"The program has encountered the following FATAL ERROR :\nProblem in opening file %s\n",File);
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
void strupr(char *chaine)
{
char *pc;

for ( pc = chaine ; *pc != 0 ; ++pc )
    *pc = toupper(*pc);
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
/* Reads a line from an .ALI file and discards comments  */
/******************************************************************************/
char read_line(FILE *in, int type, int *nb_line)
{
    char *pc;

    while ( fgets(buffer,LENBUFFER - 1,in) != NULL ) {
        ++*nb_line;
        if ( (pc = strchr(buffer,'\n')) == NULL ) {
            sprintf(buffer,"Line %d is too long (MAXIMUM = %d)",*nb_line,LENSEQ - 1);
            return(-1);
        }
        if ( *buffer != '#' ) {
            if ( (type == 0) && (*buffer != '>') ) {
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
if ( nb_spec_ali == 0 )
	{
	if ( strpbrk(buffer,"EFILPQ") != NULL )
		type_seq_ali = 1;
	else
		type_seq_ali = 0;
	}

if ( type_seq_ali == 1 )
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
    for ( nb_line = 2 , nb_spec_ali = 0 , nb_nuc_tot = 0 ;
          (drap_erreur = read_line(in_ali,0,&nb_line)) == 0 ;
          ++nb_spec_ali ){
        if ( nb_spec_ali >= MAXSPECIE ) {
            fprintf(stderr,"Maximal number of species equals %d\n",MAXSPECIE);
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
    if ( nb_spec_ali == 0 ) {
        fprintf(stderr,"This file is empty\n");
        exit(1);
    }
    lenseq_tot = MAXCARSEQ / (long) nb_spec_ali;
    if ( lenseq_tot < ((long) nb_nuc_tot + 1) ) {
        fprintf(stderr,"There are too many characters, %ld (MAX = %ld)\n",(long) nb_spec_ali * (long) (nb_nuc_tot + 1),(long) MAXCARSEQ);
        exit(1);
    }
    fclose(in_ali);


    if ( (ident_ali = (char *) malloc(nb_spec_ali * LENIDENT * sizeof(char))) == NULL ) {
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


/***********************************************************************/
void cree_dichotomie()
{
int i, niveau_parenthese;
char *pc, c1, c2;

buffer[strlen(buffer)-1] = 0;   // remove the final comma
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
/* Checks whether a distance is negative                             */
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
/* Copies an identifier into a ident                            */
/******************************************************************************/
char traite_ident(char **pb, char **pa)
{
char *pb2, *pi;

if ( (pb2 = strchr(*pb,':')) == NULL )
	{
	affich_avert("A colon following an identifier is missing in the tree",10);
	return(-1);
	}
*pb2 = 0;
strncpy(ident_arb[nb_spec],*pb + 1,LENIDENT - 1);
/*
if ( (pc = strchr(ident_arb[nb_spec],'@')) != NULL )
	{
	if ( (pi = strchr(pc + 1,'@')) != NULL )
		{
		sscanf(pi + 1,"%d",&seq_lens[nb_spec]);
// following lines to delete what comes after the second @
//		for ( ; *pi != 0 ; ++pi )
//			*pi = '_';
		}
	else
		seq_lens[nb_spec] = -1;
	}
*/
for ( pi = ident_arb[nb_spec] + strlen(ident_arb[nb_spec]) - 1 ; *pi == '_' ; --pi )
	*pi = 0;
// only the first _ is changed to a space
if ( (pi = strchr(ident_arb[nb_spec],'_')) != NULL )
	*pi = ' ';
//while ( (pi = strchr(ident_arb[nb_spec],'_')) != NULL )
//	*pi = ' ';
//if ( (pi = strstr(ident_arb[nb_spec],"ZP ")) != NULL )
//	*(pi + 2) = '_';
sprintf(*pa,"%d",nb_spec);
while ( **pa != 0 )
	++*pa;
*pb2 = ':';
*pb = pb2;
if ( ++nb_spec >= MAXSPECIE )
	{
	sprintf(buffer,"There are too many species (MAXIMUM=%d)",MAXSPECIE);
	affich_avert(buffer,10);
	return(-1);
	}
return(traite_dist(pb,pa));
}

/******************************************************************************/
/* Data file read                                  */
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
		fprintf(stderr,"The file %s has not enough lines\n",file_ARB);
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

for ( nb_spec = 0 , nb_o = 0 , nb_v = 0 , nb_f = 0 , nb_neg = 0 , pb = buffer , pa = arbre ; *pb != 0 ; ++pb )
	{
	switch(*pb)
		{
		case '(' :
			++nb_o;
			*pa++ = *pb;
			if ( *(pb + 1) != '(' )
				{
				if ( traite_ident(&pb,&pa) == -1 )
					return(-1);
				}
			break;
		case ',' :
			++nb_v;
			*pa++ = *pb;
			if ( *(pb + 1) != '(' )
				{
				if ( traite_ident(&pb,&pa) == -1 )
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

if ( (p = (P_NODE) malloc (sizeof(NODE))) == NULL )
	{
	affich_avert("Insufficient memory for creation of nodes",12);
	return(NULL);
	}
p->no_ident = -2 * MAXSPECIE;
p->dist = 0.0F;
p->size_subtree = -1.0F;
p->diff_subtrees = -1.0F;
p->nb_taxon_subtree = 0;
p->nb_taxon_no_subtree = 0;
p->div_taxon = 0;
p->pere = NULL;
p->gauche = NULL;
p->droite = NULL;
return(p);
}

/******************************************************************************/
/* Node creation                                            */
/******************************************************************************/
P_NODE creer_node(char **pa, P_NODE Pere)
{
    P_NODE p;
    char *pa2, *pa3, c;

    if ( (p = malloc_node()) == NULL )
        return(NULL);

    if ( **pa == '(' ) {
/* Fills in the left column */
        ++*pa;
        if ( (p->gauche = creer_node(pa,p)) == NULL )
            return(NULL);
        if ( **pa != ',' ) {
            affich_avert("There is a node where there is only one son",12);
            return(NULL);
        }

/* Fills in the right column  */
        ++*pa;
        if ( (p->droite = creer_node(pa,p)) == NULL )
            return(NULL);
        if ( **pa != ')' ) {
            affich_avert("There is a node without )",12);
            return(NULL);
        }

/* Fills in the node  */
        *pa += 2;
        for ( pa3 = *pa ; *pa3 != 0 ; ++pa3 )
            if ( (*pa3 == ',') || (*pa3 == ')') ) {
                c = *pa3;
                *pa3 = 0;
                break;
            }
        p->dist = (float) atof(*pa);
        p->no_ident = -2 -1 * ++nb_node;
        *pa3 = c;
        *pa = pa3;
    }
    else {
        *(pa2 = strchr(*pa,':')) = 0;
        p->no_ident = atoi(*pa);
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
/* Tree construction                                       */
/******************************************************************************/
P_NODE compute_tree()
{
    P_NODE p1, p2, p3, p4;
    char *pa;

	nb_node = 0;
    percent = 1.0F;
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
        p4->no_ident = -2;

        proot->pere = NULL;
        proot->gauche = p4;
        proot->droite = p3;
        proot->no_ident = -1;

        p4->dist = (p3->dist) * percent;
        p3->dist = p3->dist * (1.0 - percent);
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
/* Calculates the size of each subtree                       */
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
/* Removes a tree by freeing allocated memory                        */
/******************************************************************************/
void verif_tree(P_NODE pnode)
{
fprintf(stdout,"\nNoeud=%3d Dist=%.2f Size=%.2f",pnode->no_ident,pnode->dist,pnode->size_subtree); fflush(stdout);
if ( pnode->gauche == NULL )
	{
	fprintf(stdout," Species=%s",ident_arb[pnode->no_ident]); fflush(stdout);
	}

if ( pnode->pere != NULL )
	{
	if ( pnode->gauche != NULL )
		{
		fprintf(stdout," Pere=%3d Gauche=%3d Droite=%3d",(pnode->pere)->no_ident,(pnode->gauche)->no_ident,(pnode->droite)->no_ident); fflush(stdout);
		}
	else
		{
		fprintf(stdout," Pere=%3d",(pnode->pere)->no_ident); fflush(stdout);
		}
	}
else
	{
	if ( pnode->gauche != NULL )
		{
		fprintf(stdout," Gauche=%3d Droite=%3d",(pnode->gauche)->no_ident,(pnode->droite)->no_ident); fflush(stdout);
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
	for ( pi = ident_arb[pnode->no_ident] , pb = buffer , i = 0 ; i < LENIDENT - 1 ; ++pi , ++pb , ++i )
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
/* Saves the trees                                                   */
/******************************************************************************/
void save_tree_arb(char *name_DIR)
{
FILE *out;
int i;

sprintf(buffer,"%s/%s",name_DIR,file_ARB);
out = ffopen(buffer,"w");
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
/* Roots the tree at the left of the active node                                 */
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
//printf("\np0=%d p1=%d p2=%d\n",p0->no_ident,p1->no_ident,p2->no_ident);
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
        (pr->gauche)->dist = f * (1.0F - percent);
        (pr->droite)->dist = f * percent;
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
void lookfor_internal_BL(P_NODE pnode, int *n)
{
if ( (pnode->droite != NULL) && (pnode->gauche != NULL) )
	{
	internal_BL[(*n)++] = pnode->dist;
	lookfor_internal_BL(pnode->droite,n);
	lookfor_internal_BL(pnode->gauche,n);
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
//fprintf(stdout,"Noeud=%3d Dist=%.2f Size=%.2f pcour=%d\n",pnode->no_ident,pnode->dist,pnode->size_subtree,(int) pcour); fflush(stdout);
if ( pcour != NULL )
	return;
if ( (pnode->droite != NULL) && (pnode->gauche != NULL) )
	{
	lookfor_node(pnode->droite,i);
	lookfor_node(pnode->gauche,i);
	}
if ( pnode->no_ident == i )
	pcour = pnode;
}

/******************************************************************************/
/* Roots at the left of the node of the species with no_ident 0     */
/******************************************************************************/
void root_left_species0()
{
pcour = NULL;
lookfor_node(proot,0);
//fprintf(stdout,"\nNoeud=%3d Dist=%.2f Size=%.2f",pcour->no_ident,pcour->dist,pcour->size_subtree); fflush(stdout);
root_left();
}

/******************************************************************************/
/* Finds the node for midpoint rooting                       */
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
//fprintf(stdout,"Node %d: MADDIV=%d #taxa1=%d #taxa2=%d (max=%d) (diff_size= %f/%f =%f)\n",pnode->no_ident,max_div,nb_taxa1,nb_taxa2,nb_species_taxon,pnode->diff_subtrees,min_diff_subtrees,100.0*pnode->diff_subtrees/min_diff_subtrees);
	}
//fprintf(stdout,"Node=%d (b=%f ecart1-2-3=%f %f %f (min=%f)\n",pnode->no_ident,pnode->dist,f1,f2,f3,min_diff_subtrees);
}

/******************************************************************************/
/* Generates the list of all species in a subtree          */
/****************************************************************************/
void lookfor_species_subtree(P_NODE pnode)
{
if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	lookfor_species_subtree(pnode->gauche);
	lookfor_species_subtree(pnode->droite);
	}
else
	no_spec_subtree[nb_spec_subtree++] = pnode->no_ident;
}

/******************************************************************************/
/* Searches in a subtree if a species is present in the no_spec_subtree list */
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
	strcpy(buffer,ident_arb[pnode->no_ident]);
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
void generate_liste_species(P_NODE pnode)
{
int i, l;

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	generate_liste_species(pnode->gauche);
	generate_liste_species(pnode->droite);
	}
else
	{
	strcpy(buffer,ident_arb[pnode->no_ident]);
	if ( strchr(buffer,'@') != NULL )
		*strchr(buffer,'@') = 0;
	l = strlen(buffer);
	for ( i = 0 ; i < nb_species_taxon ; ++i )
//		if ( strncmp(buffer,species[i],MAX(l,strlen(species[i]))) == 0 )
		if ( strncmp(buffer,species[i],l) == 0 )
			break;
	if ( i == nb_species_taxon )
		strcpy(species[nb_species_taxon++],buffer);
	}
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
	strcpy(buffer,ident_arb[pnode->no_ident]);
	if ( strchr(buffer,'@') != NULL )
		*strchr(buffer,'@') = 0;
	l = strlen(buffer);
	for ( i = 0 ; i < nb_species_taxon ; ++i )
		if ( strncmp(buffer,species[i],l) == 0 )
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
	strcpy(buffer,ident_arb[pnode->no_ident]);
	if ( strchr(buffer,'@') != NULL )
		*strchr(buffer,'@') = 0;
	l = strlen(buffer);
	for ( i = 0 ; i < nb_species_taxon ; ++i )
		if ( strncmp(buffer,species[i],l) == 0 )
			break;
	if ( drap_species[i] == -1 )
		{
		drap_species[i] = 0;
		++nb_div_no_subtree;
		}
	}
}

/****************************************************************************/
void find_best_split_div_taxon_int(P_NODE pnode, FILE *out_stat)
{
int i;

for ( i = 0 ; i < nb_species_taxon ; ++i )
	drap_species[i] = -1;
nb_div_subtree = 0;
estimate_div_subtree(pnode);

for ( i = 0 ; i < nb_species_taxon ; ++i )
	drap_species[i] = -1;
nb_div_no_subtree = 0;
estimate_div_no_subtree(proot,pnode);

pnode->nb_taxon_subtree = nb_div_subtree;
pnode->nb_taxon_no_subtree = nb_div_no_subtree;

fprintf(out_stat,"Node%d\t%d\t%d\t%d\t%f\t%d\n",pnode->no_ident,nb_div_subtree,nb_div_no_subtree,abs(nb_div_subtree - nb_div_no_subtree),pnode->dist,nb_div_subtree + nb_div_no_subtree);

//fprintf(stdout,"tot=%d Node %d (%f): #subtree=%d #nosubtree=%d\n",nb_div_subtree+nb_div_no_subtree,pnode->no_ident,pnode->dist,nb_div_subtree,nb_div_no_subtree);

if ( max_div <= (nb_div_subtree + nb_div_no_subtree) )
	{
//fprintf(stderr,"Node %d (%f): ",pnode->no_ident,pnode->dist);
	if ( max_div == (nb_div_subtree + nb_div_no_subtree) )
		{
//fprintf(stderr,"trop court\n");
		if ( pnode->dist < max_internal_BL )
			return;
		}
	pcour = pnode;
	max_div = nb_div_subtree + nb_div_no_subtree;
	max_internal_BL = pnode->dist;
	nb_taxa1 = nb_div_subtree;
	nb_taxa2 = nb_div_no_subtree;
//fprintf(stdout,"MAXDIV=%d MAXINTERNAL=%f #taxa1=%d #taxa2=%d (max=%d) (diff_size= %f/%f =%f)\n",max_div,max_internal_BL,nb_taxa1,nb_taxa2,nb_species_taxon,pnode->diff_subtrees,min_diff_subtrees,100.0*pnode->diff_subtrees/min_diff_subtrees);
	}
}

/****************************************************************************/
void find_best_split_div_taxon(P_NODE pnode, FILE *out_stat)
{
find_best_split_div_taxon_int(pnode,out_stat);

if ( (pnode->gauche != NULL) && (pnode->droite != NULL) )
	{
	find_best_split_div_taxon(pnode->gauche,out_stat);
	find_best_split_div_taxon(pnode->droite,out_stat);
	}
}

/****************************************************************************/
int compute_nb_shared_taxa_2subtrees(void)
{
int i, j, nb_shared_taxa, nb_spec_subtree_gauche, nb_duplicates;
char ident_gauche[MAXSPECIE][LENIDENT], drap_shared[MAXSPECIE];

nb_spec_subtree = 0;
lookfor_species_subtree(proot->gauche);
nb_duplicates = 0;
for ( i = 0 ; i < nb_spec_subtree ; ++i )
	{
// To have only unique species
	for ( j = 0 ; j < (i - nb_duplicates) ; ++j )
		if ( strncmp(ident_gauche[j],ident_arb[no_spec_subtree[i]],strlen(ident_gauche[j])) == 0 )
			break;
	if ( j == (i - nb_duplicates) )
		{
		strcpy(ident_gauche[i - nb_duplicates],ident_arb[no_spec_subtree[i]]);
		if ( strchr(ident_gauche[i - nb_duplicates],'@') != NULL )
			*strchr(ident_gauche[i - nb_duplicates],'@') = 0;
		}
	else
		{
//		fprintf(stdout,"-----> %s is present multiple times in the left subtree\n",ident_arb[no_spec_subtree[i]]);
		++nb_duplicates;
		}
//fprintf(stdout,"ident_gauche[%d]=%s, ident_arb[no_spec_subtree[%d]]=%s\n",j,ident_gauche[j],i,ident_arb[no_spec_subtree[i]]);
	}
nb_spec_subtree_gauche = nb_spec_subtree - nb_duplicates;
//for ( j = 0 ; j < nb_spec_subtree_gauche ; ++j )
//	fprintf(stdout,"GAUCHE: %s\n",ident_gauche[j]);

nb_spec_subtree = 0;
lookfor_species_subtree(proot->droite);

for ( i = 0 ; i < nb_spec_subtree_gauche ; ++i )
	drap_shared[i] = -1;
for ( i = 0 ; i < nb_spec_subtree ; ++i )
	for ( j = 0 ; j < nb_spec_subtree_gauche ; ++j )
		if ( strncmp(ident_gauche[j],ident_arb[no_spec_subtree[i]],strlen(ident_gauche[j])) == 0 )
			{
			drap_shared[j] = 0;
//			++nb_shared_taxa;
			break;
			}
for ( i = 0 , nb_shared_taxa = 0 ; i < nb_spec_subtree_gauche ; ++i )
	if ( drap_shared[i] == 0 )
		++nb_shared_taxa;

return(nb_shared_taxa);
}

/****************************************************************************/
void write_split_files(char *name_DIR)
{
int i, j, nb0, nb1, nb2, nb_added_A, nb_added_B;
FILE *out1, *out2, *out_db, *out_family;
float best1, best2;
char *ps;

// to identify in which file .ali the sequences should be placed
nb0 = nb1 = nb2 = 0;
nb_spec_subtree = 0;
lookfor_species_subtree(proot->gauche);
for ( i = 0 ; i < nb_spec_ali ; ++i )
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
for ( i = 0 ; i < nb_spec_ali ; ++i )
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
if ( (nb1 == 0) || (nb2 == 0) )
	{
	fprintf(stderr,"There are not a sufficient number of species in common between %s and %s\n",file_ALI,file_ARB);
	exit(1);
	}

sprintf(file_SPLIT,"%s/%s",name_DIR,file_ALI);
*strchr(file_SPLIT,'.') = 0;
sprintf(buffer,"%s-A.ali",file_SPLIT);
out1 = ffopen(buffer,"w");
sprintf(buffer,"%s-B.ali",file_SPLIT);
out2 = ffopen(buffer,"w");
sprintf(buffer,"%s.ali",file_SPLIT);
out_family = ffopen(buffer,"w");
fprintf(out1,"#%s\n#%s\n",molecule,domaine);
fprintf(out2,"#%s\n#%s\n",molecule,domaine);
fprintf(out_family,"#%s\n#%s\n",molecule,domaine);
out_db = ffopen("db_tmp","w");
for ( i = 0 ; i < nb_spec_ali ; ++i )
	if ( drap_species[i] == 1 )
		{
		fprintf(out_family,">A-%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
		fprintf(out1,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
		fprintf(out_db, ">db1|%d\n", i);
		for ( ps = seq + i * lenseq_tot ; *ps != 0 ; ++ps )
			if ( strchr("*-$ X?",*ps) == NULL )
				fprintf (out_db, "%c",*ps);
		fprintf(out_db, "\n");
		}
	else
		if ( drap_species[i] == 2 )
			{
			fprintf(out_family,">B-%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
			fprintf(out2,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
			fprintf(out_db, ">db2|%d\n", i);
			for ( ps = seq + i * lenseq_tot ; *ps != 0 ; ++ps )
				if ( strchr("*-$ X?",*ps) == NULL )
					fprintf(out_db, "%c",*ps);
			fprintf(out_db, "\n");
			}
fclose(out1);
fclose(out2);
fclose(out_db);

//putting the SEQUENCES OF THE .ALI that are not present in the .arbre
nb_added_A = nb_added_B = 0;
for ( i = 0 ; i < nb_spec_ali ; ++i )
	if ( drap_species[i] == 0 )
		{
		if ( nb0++ == 0 )
			{
			if ( type_seq_ali == 1 )
				system("formatdb -i db_tmp -p T -o F");
			else
				system("formatdb -i db_tmp -p F -o F");
			}
		out1 = ffopen("seed.tmp","w");
		fprintf(out1, ">seed%d\n", i);
		for ( ps = seq + i * lenseq_tot ; *ps != 0 ; ++ps )
			if ( strchr("*-$ X?",*ps) == NULL )
				fprintf (out1, "%c",*ps);
		fprintf (out1, "\n");
		fclose(out1);
		if ( type_seq_ali == 1 )
			strcpy(buffer,"blastp -db db_tmp -query seed.tmp -out blastp.tmp -seg no");
		else
			strcpy(buffer,"blastn -db db_tmp -query seed.tmp -out blastp.tmp");
		if ( system(buffer) != 0 )
			fprintf(stderr,"Error in blasting %s for file %s\n",ident_ali + i * LENIDENT,file_ALI);
		else
			{
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
				fprintf(out_family,">A-%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
				sprintf(buffer,"%s-A.ali",file_SPLIT);
				out1 = ffopen(buffer,"a");
				fprintf(out1,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
				fclose(out1);
				++nb1;
				++nb_added_A;
//fprintf(stdout,"sequence absent from .arb and put in -1.ali:%s\n",ident_ali + i * LENIDENT);
				}
			else
				{
				fprintf(out_family,">B-%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
				sprintf(buffer,"%s-B.ali",file_SPLIT);
				out2 = ffopen(buffer,"a");
				fprintf(out2,">%s\n%s\n",ident_ali + i * LENIDENT,seq + (long) i * lenseq_tot);
				fclose(out2);
				++nb2;
				++nb_added_B;
//fprintf(stdout,"sequence absent from .arb and put in -2.ali:%s\n",ident_ali + i * LENIDENT);
				}
			}
		}
fclose(out_family);

if ( nb1 < nb2 )
	{
	sprintf(buffer,"mv %s-A.ali detect-problems-arb.TMP",file_SPLIT);
	system(buffer);
	sprintf(buffer,"mv %s-B.ali %s-A.ali",file_SPLIT,file_SPLIT);
	system(buffer);
	sprintf(buffer,"mv detect-problems-arb.TMP %s-B.ali",file_SPLIT);
	system(buffer);
	i = nb1;
	nb1 = nb2;
	nb2 = i;
	fprintf(stdout," (%d and %d sequences added by blast to families B and A, respectively)\n",nb_added_B,nb_added_A);
	}
else
	fprintf(stdout," (%d and %d sequences added by blast to families A and B, respectively)\n",nb_added_A,nb_added_B);
}

/****************************************************************************/
void usage(void)
{
fprintf(stderr,"USAGE: root-max-div-taxon");
fprintf(stderr,"\n\tarb=<name of file .arb>");
fprintf(stderr,"\n\tali=<name of file .ali>");
fprintf(stderr,"\n\tcriterion0=<minimum number of species in clade 1,minimum number of species in clade 2,percent of internal branches longer than the one used to split the tree,minimum number of species shared by the two subtrees> (e.g. criterion1=20,15,10,0)");
fprintf(stderr,"\n\tcriterion1=<minimum number of species in clade 1,minimum number of species in clade 2,percent of internal branches longer than the one used to split the tree,minimum number of species shared by the two subtrees> (e.g. criterion1=20,15,10,15)");
fprintf(stderr,"\n\t...\n\n");
exit(1);
}

/****************************************************************************/
int main(int argc, char **argv)
{
int i, j, n1, n2, nb_shared_taxa;
FILE *out_stat;
char *pc, name_dir[LENFILENAME], drap_SPLIT;
float f;

if ( argc < 4 )
	usage();

nb_criteria = 0;
for( i = 0 ; i < MAXCRITERIA ; ++i )
	{
	min_nb_species_group1[i] = 0;
	min_nb_species_group2[i] = 0;
	percent_longer_internal_branch[i] = 0.0;
	min_shared_taxa[i] = 0;
	}

for( i = 1 ; i < argc ; ++i )
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
			if ( (strncmp(argv[i],"criterion",9) == 0) && ((pc = strchr(argv[i] + 9,'=')) != NULL ) )
				{
				sscanf(argv[i],"criterion%d=",&j);
				sscanf(pc,"=%d,%d,%f,%d",&min_nb_species_group1[j],&min_nb_species_group2[j],&percent_longer_internal_branch[j],&min_shared_taxa[j]);
				++nb_criteria;
				}
			else
				{
				fprintf(stderr,"This argument is not correctly formatted: %s\n\n",argv[i]);
				usage();
				}
if ( nb_criteria == 0 )
	{
	fprintf(stderr,"You should provide at least one criterion\n\n");
	exit(1);
	}
for ( i = 0 ; i < nb_criteria ; ++i )
	{
	if ( (min_nb_species_group1[i] == 0) || (min_nb_species_group2[i] == 0) || (percent_longer_internal_branch[i] == 0.0) || (min_shared_taxa[i] < 0) )
		{
		fprintf(stderr,"Criterion %d is not valid: %d,%d,%f,%d\n\n",i,min_nb_species_group1[i],min_nb_species_group2[i],percent_longer_internal_branch[i],min_shared_taxa[i]);
		exit(1);
		}
	}

for ( i = 0 ; i < nb_criteria ; ++i )
	{
	sprintf(buffer,"if ! [ -d split-more%d-more%d-longer%.0f-shared%d ]; then mkdir split-more%d-more%d-longer%.0f-shared%d ; fi",min_nb_species_group1[i],min_nb_species_group2[i],percent_longer_internal_branch[i],min_shared_taxa[i],min_nb_species_group1[i],min_nb_species_group2[i],percent_longer_internal_branch[i],min_shared_taxa[i]);
	system(buffer);
	system("if ! [ -d unsplit ]; then mkdir unsplit ; fi");
	system("if ! [ -d noparalog ]; then mkdir noparalog ; fi");
	}

read_file_ALI();
if ( read_File_ARB() == -1 )
	exit(-1);

if ( (proot = compute_tree()) == NULL )
	exit(-1);
// Roots on species 0, with one of its two branches at 0, to ensure the midlength search
root_left_species0();
//verif_tree(proot);
//sprintf(file_ARB,"root0.arb");
//save_tree_arb();
compute_size_subtree(proot);

// Searches for midpoint rooting and stores the size differences between the two subtrees
min_diff_subtrees = proot->size_subtree;
pcour = proot;
lookfor_midlength(proot);
percent = 0.0;
root_left();

//Finds the node that maximizes taxonomic diversity on each side of its branch
nb_species_taxon = 0;
max_internal_BL = 0.0;
generate_liste_species(proot);

strcpy(buffer,file_ARB);
*strchr(buffer,'.') = 0;
strcat(buffer,".log");
out_log = ffopen(buffer,"w");
if ( nb_spec == nb_species_taxon )
	{
	fprintf(out_log,"There are no duplicated sequences in file %s\n",file_ALI);
	fprintf(stdout,"There are no duplicated sequences in file %s\n",file_ALI);
	sprintf(buffer,"cp %s noparalog",file_ALI);
	system(buffer);
	sprintf(buffer,"cp %s noparalog",file_ARB);
	system(buffer);
	fclose(out_log);
	exit(0);
	}

fprintf(out_log,"List of the %d unique taxa (for %d sequences):\n",nb_species_taxon,nb_spec);
for ( i = 0 ; i < nb_species_taxon ; ++i )
	fprintf(out_log,"%s\n",species[i]);

nb_internal_branch = 0;
lookfor_internal_BL(proot,&nb_internal_branch);
qsort(internal_BL,(size_t) nb_internal_branch,sizeof(float),(const void *) comparer_float);
--nb_internal_branch;		// To correct the fact that we have a rooted tree

pcour = proot;
max_div = 0;
nb_taxa1 = nb_taxa2 = 0;
strcpy(buffer,file_ARB);
*strchr(buffer,'.') = 0;
strcat(buffer,".stat");
out_stat = ffopen(buffer,"w");
fprintf(out_stat,"Node\t#species1\t#species2\tD#species1,2\tBL\tS#species1,2\n");
find_best_split_div_taxon(proot->gauche,out_stat);
find_best_split_div_taxon(proot->droite,out_stat);
fclose(out_stat);

drap_SPLIT = -1;
if ( nb_taxa1 >= nb_taxa2 )
	{
	n1 = nb_taxa1;
	n2 = nb_taxa2;
	}
else
	{
	n1 = nb_taxa2;
	n2 = nb_taxa1;
	}
fprintf(out_log,"Branch lengths of the tree sorted by decreasing value:\n");
for ( i = 0 ; i < nb_internal_branch ; ++i )
	fprintf(out_log,"%.4f ",internal_BL[i]);
fprintf(out_log,"\n");
fprintf(out_log,"\tnb_taxa(subtree 1)\tnb_taxa(subtree 2)\tBL of the best internal branch\n");
fprintf(out_log,"BEST BRANCH\t%d\t%d\t%.4f\n",n1,n2,max_internal_BL);
for ( i = 0 ; i < nb_criteria ; ++i )
	{
	j = (int) (percent_longer_internal_branch[i] * (float) nb_internal_branch / 100.0);
//fprintf(stdout,"i=%d j=%d internal_BL[j]=%f\n",i,j,internal_BL[j]);
	fprintf(out_log,"criterion%d\t%d\t%d\t%.4f\t",i,min_nb_species_group1[i],min_nb_species_group2[i],internal_BL[j]);
	if ( (n1 >= min_nb_species_group1[i]) && (n2 >= min_nb_species_group2[i]) && (max_internal_BL >= internal_BL[j]) )
		{
		sprintf(name_dir,"split-more%d-more%d-longer%.0f-shared%d",min_nb_species_group1[i],min_nb_species_group2[i],percent_longer_internal_branch[i],min_shared_taxa[i]);
		root_left();
		percent = 0.5;
		f = (proot->gauche)->dist + (proot->droite)->dist;
		(proot->gauche)->dist = (1.0 - percent) * f;
		(proot->droite)->dist = percent * f;
		nb_shared_taxa = compute_nb_shared_taxa_2subtrees();
		fprintf(out_log,"nb taxa shared by the 2 subtrees\t%d\n",nb_shared_taxa);
		if ( nb_shared_taxa >= min_shared_taxa[i] ) 
			{
			save_tree_arb(name_dir);
//fprintf(stdout,"File %s has been split according to criterion %d: (%d-%d species) Node %d: MADDIV= %d #taxa1= %d #taxa2= %d (max= %d) (diff_size= %f/%f %.0f) (internal= %f /%f)\n",file_ALI,i,n1,n2,pcour->no_ident,max_div,pcour->nb_taxon_subtree,pcour->nb_taxon_no_subtree,nb_species_taxon,pcour->diff_subtrees,min_diff_subtrees,100.0*pcour->diff_subtrees/min_diff_subtrees,pcour->dist,internal_BL[(int) (percent_longer_internal_branch[0] * (float) nb_internal_branch)]);
fprintf(stdout,"File %s has been split according to criterion %d: (%d>=%d, %d>=%d, %.4f>=%.4f %d>=%d)",file_ALI,i,n1,min_nb_species_group1[i],n2,min_nb_species_group2[i],max_internal_BL,internal_BL[j],nb_shared_taxa,min_shared_taxa[i]);
//fprintf(stdout,"File %s has been SPLITTED (%d-%d species) Node %d: MADDIV= %d #taxa1= %d #taxa2= %d (max= %d) (diff_size= %f/%f %.0f) (internal= %f /%f)\n",file_ALI,nb1,nb2,pcour->no_ident,max_div,pcour->nb_taxon_subtree,pcour->nb_taxon_no_subtree,nb_species_taxon,pcour->diff_subtrees,min_diff_subtrees,100.0*pcour->diff_subtrees/min_diff_subtrees,pcour->dist,internal_BL[(int) (percent_longer_internal_branch * (float) nb_internal_branch)]);
			write_split_files(name_dir);
			drap_SPLIT = 0;
			fprintf(out_log,"OK\n");
			}
		else
			fprintf(out_log,"NO because of too few shared taxa\n");
		}
	else
		fprintf(out_log,"NO\n");
	}
if ( drap_SPLIT == -1 )
	{
	percent = 0.5;
	f = (proot->gauche)->dist + (proot->droite)->dist;
	(proot->gauche)->dist = (1.0 - percent) * f;
	(proot->droite)->dist = percent * f;
	save_tree_arb("unsplit");
	sprintf(buffer,"cp %s unsplit",file_ALI);
	system(buffer);
	fprintf(stdout,"File %s has not been split\n",file_ALI);
	}
fclose(out_log);
system("rm -f blastp.tmp db_tmp db_tmp.p?? formatdb.log seed.tmp");
return(0);
}
