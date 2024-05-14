/* The PITESOFT computes the Primary Indirect Topographic Effect in Stokes-Helmert method.
 * This program is distributed under the terms of GNU License published by the author.
 * Hopefully, this program will be useful without warranty. 
 *
 * Author    : R. Alpay ABBAK
 * Address   : Konya Technical University, Geomatics Engineering Dept, Campus, Konya, TURKEY
 * E-mail    : raabbak@ktun.edu.tr
 *
 * Please cite the following paper for the program:
 * Oztop, A, Abbak RA, Ustun A (2021). The effect of crustal density variation on PITE 
 * (Primary Indirect Topography Effect), Geomatik, 9 (1), 97-105. 
 *
 * Compilation of the program on Linux:
 * g++ PITESOFT.cpp matris.cpp -o PITESOFT
 *
 * Execution of the program with the sample data: 
 * ./PITESOFT -Eelevation.xyz -Ddensity.xyz -R45.005/47/2.005/4 -I0.01/0.01
 * Created : 25.01.2022				v1.0  
 * Updated : 14.05.2024				v2.0  */
#include <unistd.h>				// Standard option library 
#include "matris.h"				// User-defined matrix library
#define G 6.6742E-11				// Newtonian Earth's attraction constant
#define R 6371000.0				// Earth mean radius (m)
#define pi 3.141592653589793			// constant pi
#define rho 57.2957795130823			// 180degree/pi
#define rho0 2670.0 				// mean crust density 
#define mx 800					// Maximum dimension of the data area 
#define OPTIONS "D:E:P:R:I:H"			// Parameters and Options
void HELP();
int main(int argc, char *argv[])
{//*** D E F I N I T I O N S **********************************************************************/
	FILE *elevation=NULL;			// elevation data file
	FILE   *density=NULL;			// density data file
	char 	      option;			// the whole option line
    	const char *slash="/";			// discriminant symbol
    	int           i=0;			// array element
    	int           j=0;			// array element
	int           x=0;			// array element
    	int           y=0;			// array element
    	int          in=0;			// array element of computation point
    	int          jn=0;			// array element of computation point
	double 	 MinPhi=45.005;			// minimum latitude of target area
	double 	 MaxPhi=47.00;			// maximum latitude of target area
	double 	 MinLam=02.005;			// minimum longitude of target area
	double 	 MaxLam=04.00;			// maximum longitude of target area
	double 	 PhiInt=0.01;			// grid size in latitude direction
	double 	 LamInt=0.01;			// grid size in longitude direction
	double MinDatPhi=0.0;			// minimum latitude of data area
	double MinDatLam=0.0;			// minimum longitude of data area
	double      phi=0.0;			// latitude of computational point [deg] 
    	double      lam=0.0;			// longitude of computational point [deg]
    	double   total1=0.0;			// total value
    	double   total2=0.0;			// total value
    	double   total3=0.0;			// total value
    	double   total4=0.0;			// total value
    	double     bel0=0.0;			// total value
    	double     bel1=0.0;			// total value
    	double     bel2=0.0;			// total value
	double    gamma=0.0;			// normal gravity on ellipsoid
	double 	      l=0.0;
	double psi0=1.5/rho;			// integration capsize
	double       r1=0.0;			// HP+R
	double       r2=0.0;			// HQ+R
	double       cp=0.0;			// cos(psi) 
	matris 	   lati(mx);			// latitude of running point
	matris     loni(mx);			// longitude of running point
	matris     H(mx,mx);			// topographic heights 
	matris     D(mx,mx);			// topographic densities 
	matris    NH(mx,mx);			// topographic correction on geoid
	matris   NH0(mx,mx);			// topographic correction on geoid
	matris   NH1(mx,mx);			// topographic correction on geoid
	matris   NH2(mx,mx);			// topographic correction on geoid
	matris   NH3(mx,mx);			// topographic correction on geoid
	matris   NH4(mx,mx);			// topographic correction on geoid
	matris   psi(mx,mx);			// computational psi between P and Q
/***** O P T I O N S   A N A L Y S I S ***********************************************************/
	if(argc<2) 
	{
		HELP();    exit(EXIT_FAILURE); 
	}
	while((option=getopt(argc,argv,OPTIONS))!=-1)
        switch(option)
        {
		case 'D':
			density=fopen(optarg,"r");	
                	break;
		case 'E':
			elevation=fopen(optarg,"r");	
                	break;
		case 'P':
			psi0=atof(optarg)/rho;
			break;
		case 'R':
			MinPhi=atof(strtok(optarg,slash));
			MaxPhi=atof(strtok(NULL,slash));
			MinLam=atof(strtok(NULL,slash));
			MaxLam=atof(strtok(NULL,slash));
                	break;
		case 'I':
                	PhiInt=atof(strtok(optarg,slash));
                	LamInt=atof(strtok(NULL,slash));
                	break;
		case 'H':
                	HELP();
                	exit(EXIT_SUCCESS);
		default:
                	HELP();
                	exit(EXIT_FAILURE);
        }
//***** R E A D I N G   F I L E S   &  D E F I N I N G   L I M I T S ******************************/
	if(elevation==NULL)
    	{
		printf("\nElevation file can not be opened!!!\n\n");
        	exit(EXIT_FAILURE);
    	}
	MinDatPhi=MinPhi-3.0;
	MinDatLam=MinLam-3.0;
	while(!feof(elevation))
    	{
		fscanf(elevation,"%lf%lf",&phi,&lam);
		i=round((phi-MinDatPhi)/PhiInt); 	j=round((lam-MinDatLam)/LamInt);
		lati(i)=phi/rho; 			loni(j)=lam/rho;
		fscanf(elevation,"%lf\n",&H(i,j));
    	}
    	fclose(elevation);
	while(!feof(density))
    	{
		fscanf(density,"%lf%lf",&phi,&lam);
		i=round((phi-MinDatPhi)/PhiInt); 	j=round((lam-MinDatLam)/LamInt);
		fscanf(density,"%lf\n",&D(i,j));
    	}
    	fclose(density);
	int     framePhi=round(psi0*rho/PhiInt);	// vertical limit of compartment
    	int     frameLam=round(acos((cos(psi0)-pow(sin(MaxPhi/rho),2))/(pow(cos(MaxPhi/rho),2)))*rho/LamInt);//horizontal limit	
	double    const1=4.0*pi*G*rho0;			// grid size of the block 
	double    const2=G*rho0*PhiInt*LamInt/rho/rho; 	// coefficient for the second part of NH
	double    const3=G*PhiInt*LamInt/rho/rho; 	// coefficient for the second part of NH
/***** C O M P U T E   S E C T I O N ************************************************************/
	for(phi=MinPhi;phi<=MaxPhi;phi=phi+PhiInt)
   	{
		lam=MinLam;
		in=round((phi-MinDatPhi)/PhiInt);
		jn=round((lam-MinDatLam)/LamInt);
		for(i=in-framePhi;i<=in+framePhi;i++)
		{
			x=i-in+framePhi; 			// position on the compartment
			for(j=jn;j<=jn+frameLam;j++)
			{
				y=abs(j-jn); 			// position on the compartment
				if(i==in && j==jn) psi(x,y)=0.00001;
				else psi(x,y)=acos(sin(lati(in))*sin(lati(i))+cos(lati(in))*cos(lati(i))*cos(loni(j)-loni(jn)));
			}
		}
	        gamma=9.7803267715*(1.0+0.001931851353*pow(sin(phi/rho),2))/sqrt(1.0-0.006694380023*pow(sin(phi/rho),2));		
		for(lam=MinLam;lam<=MaxLam;lam=lam+LamInt)
		{
			in=round((phi-MinDatPhi)/PhiInt); // computation point
			jn=round((lam-MinDatLam)/LamInt); // computation point
			r1=H(in,jn)+R;  		  // r_P
			for(i=in-framePhi;i<=in+framePhi;i++)
			{
				x=i-in+framePhi; 	  // position on the compartment
				for(j=jn-frameLam;j<=jn+frameLam;j++)
				{
					y=abs(j-jn); 		// position on the compartment
					if(psi(x,y)<psi0)       // 
					{
						r2=H(i,j)+R;	// r_Q
						cp=cos(psi(x,y));
						l=sqrt(R*R+r2*r2-2.0*R*r2*cp);
						bel0=(r1+3.0*r1*cp)*l+(r1*r1*(3.0*cp*cp-1.0)*log(r1-r1*cp+l));
						bel1=(r2+3.0*r1*cp)*l+(r1*r1*(3.0*cp*cp-1.0)*log(r2-r1*cp+l));
						bel2=(R+3.0*r1*cp)*l+(r1*r1*(3.0*cp*cp-1.0)*log(R-r1*cp+l));
						total1+=0.5*(bel1-bel0)*cos(lati(i));
						total2+=(pow(r2,3)-pow(r1,3))/3.0/l*cos(lati(i));
						total3+=(2670.0-D(i,j))*0.5*(bel2-bel1)*cos(lati(i));
						total4+=(2670.0-D(i,j))*(pow(r2,3)-pow(R,3))/3.0/l*cos(lati(i));
					}
				}
			}
			NH0(in,jn)=-const1/gamma*H(in,jn)*H(in,jn)*(0.5+H(in,jn)/3.0/R);
			NH1(in,jn)= const2/gamma*total1;
			NH2(in,jn)=-const2/gamma*total2;
			NH3(in,jn)= const3/gamma*total3;
			NH4(in,jn)= const3/gamma*total4;
			NH(in,jn)=NH0(in,jn)+NH1(in,jn)+NH2(in,jn)+NH3(in,jn)+NH4(in,jn);
			printf("%9.6f%9.6f%9.4f%9.4f%9.4f%9.4f%9.4f%9.4f\n",phi,lam,NH0(in,jn),NH1(in,jn),NH2(in,jn),NH3(in,jn),NH4(in,jn),NH(in,jn));
			total1=total2=total3=total4=0.0;
		}
	}
	return 0;
/***** F I N I S H *******************************************************************************/
}
void HELP()
{
	fprintf(stderr,"\n     The PITESOFT computes the Primary Indirect Topographic Effect in Stokes-Helmert method.\n\n");
	fprintf(stderr,"USAGE:\n");
	fprintf(stderr,"     PITESOFT -E[elevation<file>] -D[density<file>] -P[<value>] -R[<value>] -I[<value>] -H\n\n");
	fprintf(stderr,"PARAMETERS:\n");
	fprintf(stderr,"     elevation: mean topographic elevations which cover the data area.\n");
	fprintf(stderr,"                It includes grid based DEM data (latitude, longitude, and mean elevation).\n\n");
	fprintf(stderr,"     density:   topgraphic densities which cover the data area.\n");
	fprintf(stderr,"                It includes grid based data (latitude, longitude, and density).\n\n");
	fprintf(stderr,"OPTIONS:\n");
	fprintf(stderr,"     -P<value>  integration capsize (unit: degree).\n");
	fprintf(stderr,"                default: 1.5\n\n");
	fprintf(stderr,"     -R<value>  limits of the target area (MinLat/MaxLat/MinLon/MaxLon).\n");
	fprintf(stderr,"                default: 45.005/46.995/02.005/3.995\n\n");
	fprintf(stderr,"     -I<value>  intervals of the grid (LatInterval/LonInterval).\n");
	fprintf(stderr,"                default: 0.01/0.01 (36x36 arc-seconds)\n\n");
	fprintf(stderr,"     -H         prints this help.\n\n");
	fprintf(stderr,"Developed by:   Dr. R. Alpay Abbak\n\n");
}
