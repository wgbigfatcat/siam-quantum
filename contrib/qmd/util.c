#include <stdlib.h>
#include <stdio.h>

#include "qmd.h"
#include "util.h"

struct AtomName_t{
	int   Z;
	char  *shortName;
	char  *longName;
};

//
// Nov 4, 2009 - Theerapon Kumla
//		Fix typo from Selinium to Selenium
//		Add atomic number 55 - 103
//
#define MAX_PERIODIC_ATOM 103
static struct AtomName_t PeriodicName[MAX_PERIODIC_ATOM] =
{
              {1, "H", "HYDROGEN"},  \
              {2, "HE","HELIUM"},    \
              {3, "LI","LITHIUM"},   \
              {4, "BE","BERYLLIUM"}, \
              {5, "B", "BORON"},     \
              {6, "C", "CARBON"},    \
              {7, "N", "NITROGEN"},  \
              {8, "O", "OXYGEN"},    \
              {9, "F", "FLUORINE"},  \
              {10,"NE","NEON"},      \
              {11,"NA","SODIUM"},    \
              {12,"MG","MAGNESIUM"}, \
              {13,"AL","ALUMINUM"},  \
              {14,"SI","SILICON"},   \
              {15,"P", "PHOSPHORUS"},\
              {16,"S", "SULFUR"},    \
              {17,"CL","CHLORINE"},  \
              {18,"AR","ARGON"},     \
              {19,"K", "POTASSIUM"}, \
              {20,"CA","CALCIUM"},   \
              {21,"SC","SCANDIUM"},  \
              {22,"TI","TITANIUM"},  \
              {23,"V", "VANADIUM"},  \
              {24,"CR","CHROMIUM"},  \
              {25,"MN","MANGANESE"}, \
              {26,"FE","IRON"},      \
              {27,"CO","COBALT"},    \
              {28,"NI","NICKEL"},    \
              {29,"CU","COPPER"},    \
              {30,"ZN","ZINC"},      \
              {31,"GA","GALLIUM"},   \
              {32,"GE","GERMANIUM"}, \
              {33,"AS","ARSENIC"},   \
              {34,"SE","SELENIUM"},  \
              {35,"BR","BROMINE"},   \
              {36,"KR","KRYPTON"},   \
              {37,"RB","RUBIDIUM"},  \
              {38,"SR","STRONTIUM"}, \
              {39,"Y", "YTTRIUM"},   \
              {40,"ZR","ZIRCONIUM"}, \
              {41,"NB","NIOBIUM"},   \
              {42,"MO","MOLYBDENUM"},\
              {43,"TC","TECHNETIUM"},\
              {44,"RU","RUTHENIUM"}, \
              {45,"RH","RHODIUM"},   \
              {46,"PD","PALLADIUM"}, \
              {47,"AG","SILVER"},    \
              {48,"CD","CADMIUM"},   \
              {49,"IN","INDIUM"},    \
              {50,"SN","TIN"},       \
              {51,"SB","ANTIMONY"},  \
              {52,"TE","TELLURIUM"}, \
              {53,"I", "IODINE"},    \
              {54,"XE","XENON"},     \
              {55,"CS","CESIUM"},    \
              {56,"BA","BARIUM"},    \
              {57,"LA","LANTHANUM"}, \
              {58,"CE","CERIUM"},    \
              {59,"PR","PRASEODYMIUM"}, \
              {60,"ND","NEODYMIUM"},    \
              {61,"PM","PROMETHIUM"},   \
              {62,"SM","SAMARIUM"},     \
              {63,"EU","EUROPIUM"},     \
              {64,"GD","GADOLIUM"},     \
              {65,"TB","TERBIUM"},      \
              {66,"DY","DYSPROSIUM"},   \
              {67,"HO","HOLMIUM"},      \
              {68,"ER","ERBIUM"},       \
              {69,"TM","THULIUM"},      \
              {70,"YB","YTTERBIUM"},    \
              {71,"LU","LUTETIUM"},     \
              {72,"HF","HAFNIUM"},      \
              {73,"TA","TANTALUM"},     \
              {74,"W","TUNGSTEN"},      \
              {75,"RE","RHENIUM"},      \
              {76,"OS","OSMIUM"},       \
              {77,"IR","IRIDIUM"},      \
              {78,"PT","PLATINUM"},     \
              {79,"AU","GOLD"},         \
              {80,"HG","MERCURY"},      \
              {81,"TL","THALLIUM"},     \
              {82,"PB","LEAD"},         \
              {83,"BI","BISMUTH"},      \
              {84,"PO","POLONIUM"},     \
              {85,"AT","ASTATINE"},     \
              {86,"RN","RADON"},        \
              {87,"FR","FRANCIUM"},     \
              {88,"RA","RADIUM"},       \
              {89,"AC","ACTINIUM"},     \
              {90,"TH","THORIUM"},      \
              {91,"PA","PROTACTINIUM"}, \
              {92,"U","URANIUM"},       \
              {93,"NP","NEPTUNIUM"},    \
              {94,"PU","PLUTONIUM"},    \
              {95,"AM","AMERICIUM"},    \
              {96,"CM","CURIUM"},       \
              {97,"BK","BERKELIUM"},    \
              {98,"CF","CALIFORNIUM"},  \
              {99,"ES","EINSTEINIUM"},  \
              {100,"FM","FERMIUM"},     \
              {101,"MD","MENDELEVIUM"}, \
              {102,"NO","NOBELIUM"},    \
              {103,"LR","LAWRENCIUM"}   
};

// sym2Z : convert atom name which could be in both short 
// or long format into atomic number. Any name that does not
// match string will be Z=0
//
// Feb 26, 2008 - Teepanis Chachiyo
//     Initial implementation
//
#ifndef SYMB_SHORTNAME
#define SYMB_SHORTNAME 0
#endif
#ifndef SYMB_LONGNAME
#define SYMB_LONGNAME  1
#endif
int sym2Z(const char *sym, int type){
	char str[1024];
	int i;

	// make sym upper case
	for(i=0; i < strlen(sym); i++)
		str[i] = toupper(sym[i]);
	// null termination
	str[i] = '\0';

	for(i=0; i < MAX_PERIODIC_ATOM; i++){
		switch(type){
		case SYMB_SHORTNAME:
			if(strcmp(str,PeriodicName[i].shortName)==0)
				return PeriodicName[i].Z;
			break;
		case SYMB_LONGNAME:
			if(strcmp(str,PeriodicName[i].longName)==0)
				return PeriodicName[i].Z;
			break;
		}
	}
	return 0;
}

// Summary:    Convert from atomic number to standard atom
//             name defined in PeriodicName[..] above
//
// Parameters: Z   - atomic number
//             sym - pointer to recieve atom name
//
// Author:     Teepanis Chachiyo
//
// Created:    2005
//
// Modified:
//
//	Nov 4, 2009 - Theerapon Kumla
//		catch error if name not found
//
void Z2sym(int Z, char *sym){
	if(0 < Z && Z <= MAX_PERIODIC_ATOM)
		strcpy(sym, PeriodicName[Z-1].shortName);
	else{
		printf("Z2sym - Error : Z = %d\n",Z);
		exit(EXIT_FAILURE);
	}
}


// Summary:    Search for a set of "n" string str.
//             If one of the string given is "*", then it will
//             match anything (i.e. wildcard).
//
// Parameters: fd  - file pointer to search
//             n   - number of search strings
//             str - list of keywords
//
// Return:     EOF if the end of file
//
// Author:     Teepanis Chachiyo
//
// Created:    08/02/2005
//
// Modified:
//
int findf(FILE *fd, int n, char *str, ...){

  int  s;
  int  ret;
  char buf[1024];

  for(s = 0; (ret=fscanf(fd, "%s", buf)) == 1; ){

    if( strncmp(buf, *(&str+s), 1024) == 0 ||
        strcmp("*", *(&str+s)) == 0) s++; else s = 0;

    if( s == n ) break;
  }

  return ret;
}

// change atomic number to mass in kg
void Z2kg(struct molecule_t *ptr){                      
	int i;                               // loop index
	
	// switch case	
	for ( i=0; i < ptr->nAtom; i++){	    
		switch(ptr->Z[i]){
			case 1 : 	ptr->m[i]=1.008; 	break;
			case 2 :	ptr->m[i]=4.003; 	break;
			case 3 : 	ptr->m[i]=6.941; 	break;
			case 4 : 	ptr->m[i]=9.012; 	break;
			case 5 : 	ptr->m[i]=10.81; 	break;
			case 6 : 	ptr->m[i]=12.01; 	break;
			case 7 : 	ptr->m[i]=14.01; 	break;
			case 8 : 	ptr->m[i]=16.00; 	break;
			case 9 : 	ptr->m[i]=19.00; 	break;
			case 10 : 	ptr->m[i]=20.18; 	break;
			case 11 : 	ptr->m[i]=22.99;	break;
			case 12 : 	ptr->m[i]=24.31; 	break;
			case 13 : 	ptr->m[i]=26.98; 	break;
			case 14 : 	ptr->m[i]=28.09; 	break;
			case 15 : 	ptr->m[i]=30.97; 	break;
			case 16 : 	ptr->m[i]=32.07; 	break;
			case 17 : 	ptr->m[i]=35.45;	break;
			case 18 : 	ptr->m[i]=39.95; 	break;
			case 19 : 	ptr->m[i]=39.10; 	break;
			case 20 : 	ptr->m[i]=40.08; 	break;
			case 21 : 	ptr->m[i]=44.96; 	break;
			case 22 : 	ptr->m[i]=47.87; 	break;
			case 23 : 	ptr->m[i]=50.94; 	break;
			case 24 : 	ptr->m[i]=52.00; 	break;
			case 25 : 	ptr->m[i]=54.94; 	break;
			case 26 : 	ptr->m[i]=55.84; 	break;
			case 27 : 	ptr->m[i]=58.93; 	break;
			case 28 : 	ptr->m[i]=58.69; 	break;
			case 29 : 	ptr->m[i]=63.55; 	break;
			case 30 : 	ptr->m[i]=65.39; 	break;
			case 31 : 	ptr->m[i]=69.72; 	break;
			case 32 : 	ptr->m[i]=72.61; 	break;
			case 33 : 	ptr->m[i]=74.92; 	break;
			case 34 : 	ptr->m[i]=78.96; 	break;
			case 35 : 	ptr->m[i]=79.90; 	break;
			case 36 : 	ptr->m[i]=83.80; 	break;
			case 37 : 	ptr->m[i]=85.47; 	break;
			case 38 :	ptr->m[i]=87.62; 	break;
			case 39 : 	ptr->m[i]=88.91; 	break;
			case 40 : 	ptr->m[i]=91.22; 	break;
			case 41 : 	ptr->m[i]=92.91; 	break;
			case 42 : 	ptr->m[i]=95.94; 	break;
			case 43 : 	ptr->m[i]=99.00; 	break;
			case 44 : 	ptr->m[i]=101.07;	break;
			case 45 : 	ptr->m[i]=102.91; 	break;
			case 46 : 	ptr->m[i]=106.42; 	break;
			case 47 : 	ptr->m[i]=107.87; 	break;
			case 48 : 	ptr->m[i]=112.41; 	break;
			case 49 :	ptr->m[i]=114.82; 	break;
			case 50 : 	ptr->m[i]=118.71; 	break;
			case 51 : 	ptr->m[i]=121.76; 	break;
			case 52 : 	ptr->m[i]=127.60; 	break;
			case 53 : 	ptr->m[i]=126.90; 	break;
			case 54 : 	ptr->m[i]=131.29; 	break;
			case 55 : 	ptr->m[i]=132.91; 	break;
			case 56 : 	ptr->m[i]=137.33; 	break;
			case 57 : 	ptr->m[i]=138.91; 	break;
			case 58 : 	ptr->m[i]=140.12; 	break;
			case 59 : 	ptr->m[i]=140.91; 	break;
			case 60 : 	ptr->m[i]=144.24; 	break;
			case 61 : 	ptr->m[i]=145.00; 	break;
			case 62 : 	ptr->m[i]=150.36; 	break;
			case 63 : 	ptr->m[i]=151.96; 	break;
			case 64 : 	ptr->m[i]=157.25; 	break;
			case 65 : 	ptr->m[i]=158.93; 	break;
			case 66 : 	ptr->m[i]=162.50; 	break;
			case 67 : 	ptr->m[i]=164.93; 	break;
			case 68 : 	ptr->m[i]=167.26;	break;
			case 69 : 	ptr->m[i]=168.93; 	break;
			case 70 : 	ptr->m[i]=173.04; 	break;
			case 71 : 	ptr->m[i]=174.97; 	break;
			case 72 : 	ptr->m[i]=178.49; 	break;
			case 73 : 	ptr->m[i]=180.95; 	break;
			case 74 : 	ptr->m[i]=183.84; 	break;
			case 75 : 	ptr->m[i]=186.21; 	break;
			case 76 : 	ptr->m[i]=190.23; 	break;
			case 77 : 	ptr->m[i]=192.22; 	break;
			case 78 : 	ptr->m[i]=195.08; 	break;
			case 79 : 	ptr->m[i]=196.97; 	break;
			case 80 : 	ptr->m[i]=200.59; 	break;
			case 81 : 	ptr->m[i]=204.38; 	break;
			case 82 :	ptr->m[i]=207.2; 	break;
			case 83 : 	ptr->m[i]=208.98; 	break;
			case 84 : 	ptr->m[i]=209; 		break;
			case 85 : 	ptr->m[i]=210; 		break;
			case 86 : 	ptr->m[i]=222; 		break;
			case 87 : 	ptr->m[i]=223; 		break;
			case 88 : 	ptr->m[i]=226; 		break;
			case 89 : 	ptr->m[i]=227; 		break;
			case 90 : 	ptr->m[i]=232.04; 	break;
			case 91 : 	ptr->m[i]=231.04; 	break;
			case 92 : 	ptr->m[i]=238.03; 	break;
			case 93 : 	ptr->m[i]=237; 		break;
			case 94 : 	ptr->m[i]=244; 		break;
			case 95 : 	ptr->m[i]=243; 		break;
			case 96 : 	ptr->m[i]=247; 		break;
			case 97 : 	ptr->m[i]=247; 		break;
			case 98 : 	ptr->m[i]=251; 		break;
			case 99 : 	ptr->m[i]=252; 		break;
			case 100 : 	ptr->m[i]=257; 		break;
			case 101 : 	ptr->m[i]=258; 		break;
			case 102 : 	ptr->m[i]=259; 		break;
			case 103 : 	ptr->m[i]=262; 		break;
		// data from http://en.wikipedia.org/wiki/Atomic_weight/Table
		}
		
		// change amu to kg
		ptr->m[i]= ptr->m[i]*PROTONMASS;
	}
}
