/*
This program simulates a double pendulum using RK4 numerical approximation.
The code uses a pipe to gnuplot to create a visual animation of the trajectory
of the pendulum.
Oli Turner

Note: The GNU pipe implementation was done on Linux, so unfortunately this version of code will not produce an animation of the 
double pendulum in Windows.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define g 9.81
//Pendulum
typedef struct {
    double m1; //mass of the first bob
    double m2; //mass of the second bob
    double l1; //length of the upper pendulum
    double l2; //length of the lower pendulum
    double t1; //Angle theta1 of the upper pendulum
    double t2; //Angle theta2 of the lower pendulum
    double w1; //Angular velocity of the upper pendulum bob
    double w2; //Angular velocity of the lower pendulum bob
}Pendulum;

//The equations of motion that need to be solved.
double dtheta1(double w1){
    return w1;
}
double dtheta2(double w2){
    return w2;
}
double dw1(double t1,double t2, double w1, double w2,double m1,double m2,double l1, double l2){
    double b = ((m2*l1*w1*w1*sin(t2-t1)*cos(t2-t1)) + (m2*g*sin(t2)*cos(t2-t1)) + (m2*l2*w2*w2*sin(t2-t1)) - ((m1+m2)*g*sin(t1)))/((m1+m2)*l1 - (m2*l1*cos(t2-t1)*cos(t2-t1)));
return b;
}
double dw2(double t1,double t2, double w1, double w2,double m1, double m2, double l1, double l2){
    double k = (-m2*l2*w2*w2*sin(t2-t1)*cos(t2-t1))  +  (m1+m2)*((g*sin(t1)*cos(t2-t1)) - (l1*w1*w1*sin(t2-t1)) - (g*sin(t2)))/((m1+m2)*l2 - (m2*l2*cos(t2-t1)*cos(t2-t1)));
return k;
}

//The implemented Runge Kutta 4th Order algorithm.
int rk4(Pendulum *p, double h){
    double k1,k2,k3,k4; //t1
    double p1,p2,p3,p4; //t2
    double z1,z2,z3,z4; //w1
    double j1,j2,j3,j4; //w2

    k1 = dtheta1(p->w1);
    p1 = dtheta2(p->w2);
    z1 = dw1(p->t1,p->t2,p->w1,p->w2,p->m1,p->m2,p->l1,p->l2);
    j1 = dw2(p->t1,p->t2,p->w1,p->w2,p->m1,p->m2,p->l1,p->l2);

    k2 = dtheta1(p->w1 + (h/2)*z1);
    p2 = dtheta2(p->w2 + (h/2)*j1);
    z2 = dw1((p->t1 + (h/2)*k1),(p->t2 + (h/2)*p1),(p->w1 + (h/2)*z1),(p->w2+(h/2)*j1),p->m1,p->m2,p->l1,p->l2);
    j2 = dw2((p->t1 + (h/2)*k1),(p->t2 + (h/2)*p1),(p->w1 + (h/2)*z1),(p->w2+(h/2)*j1),p->m1,p->m2,p->l1,p->l2);

    k3 = dtheta1(p->w1 + (h/2)*z2);
    p3 = dtheta2(p->w2 + (h/2)*j2);
    z3 = dw1((p->t1 + (h/2)*k2),(p->t2 + (h/2)*p2),(p->w1 + (h/2)*z2),(p->w2+(h/2)*j2),p->m1,p->m2,p->l1,p->l2);
    j3 = dw2((p->t1 + (h/2)*k2),(p->t2 + (h/2)*p2),(p->w1 + (h/2)*z2),(p->w2+(h/2)*j2),p->m1,p->m2,p->l1,p->l2);

    k4 = dtheta1(p->w1 + (z3)*h);
    p4 = dtheta2(p->w2 + (j3)*h);
    z4 = dw1((p->t1 + (k3)*h),(p->t2 + (p3)*h),(p->w1 + (z3)*h),(p->w2+(j3)*h),p->m1,p->m2,p->l1,p->l2);
    j4 = dw2((p->t1 + (k3)*h),(p->t2 + (p3)*h),(p->w1 + (z3)*h),(p->w2+(j3)*h),p->m1,p->m2,p->l1,p->l2);

    p->t1 += h*((k1/6)+(k2/3)+(k3/3)+(k4/6));
    p->t2 += h*((p1/6)+(p2/3)+(p3/3)+(p4/6));
    p->w1 += h*((z1/6)+(z2/3)+(z3/3)+(z4/6));
    p->w2 += h*((j1/6)+(j2/3)+(j3/3)+(j4/6));

    return 1;
}

int main(){

    double temp,scalingfactor;
    double t = 0,h=0.01; //t is total time, h is the time step. (dt).
    double x1=0,y1=0,x2=0,y2=0,x0=0,y0=0; //(x0,y0) is origin. (x1,y1) is coordinates of end of upper pendulum, (x2,y2) coordinates of lower pend.
    int choice;
    int i = 0;
    char runnable;
    Pendulum * dblpend; //Initialising a double pendulum & allocating space.
    dblpend = malloc(sizeof(Pendulum));

    //Frontend user input with NO validity checking - you COULD enter unrecommened values but how ridiculous is a pendulum with 100m long arms... (goes sort of weird)
    printf("DOUBLE PENDULUM SIMULATOR\n");
    printf("Can your system run the full program? (requiring linux and commandline gnuplot) (y/n)\n");
    scanf("%c",&runnable);

    printf("Would you like to: \n(1) Set your own parameters and initial conditions\n(2) Use preset values\n");
    printf("Please enter the number of your choice:\n");
    scanf("%i",&choice);

    switch(choice){
        case(1):
            printf("Please enter the mass of the upper bob (kg): (recommended between 0.1 - 3)\n");
            scanf("%lg",&dblpend->m1);
            printf("Please enter the mass of the lower bob (kg): (recommended between 0.1 - 3)\n");
            scanf("%lg",&dblpend->m2);
            printf("Please enter the length of the upper pendulum (m): (recommended between 0.3 - 2)\n");
            scanf("%lg",&dblpend->l1);
            printf("Please enter the length of the lower pendulum (m): (recommended between 0.3 - 2)\n");
            scanf("%lg",&dblpend->l2);
            printf("Please enter the starting angle of the upper pendulum: (In degrees between 0 and 180)\n");
            scanf("%lg",&temp);
            dblpend->t1 = temp * (M_PI/180);
            printf("Please enter the starting angle of the lower pendulum: (In degrees between 0 and 180)\n");
            scanf("%lg",&temp);
            dblpend->t2 = temp * (M_PI/180);
            dblpend->w1 = 0;
            dblpend->w2 = 0;
            break;

        case(2):
            dblpend->m1 = 1.00;
            dblpend->m2 = 0.90;
            dblpend->l1 = 1.00;
            dblpend->l2 = 0.90;
            dblpend->t1 = M_PI/2;
            dblpend->t2 = 0;
            dblpend->w1 = 0;
            dblpend->w2 = 0;
            break;

        default:
            printf("Please enter an option given");
            return 1;
            break;
    }
    //20s gives me 7mb gif files!!! wow!!! thank goodness for gifsicle - optimising that down to like 500kb!!.
    printf("How long would you like the simulation to run for?: (please enter number of seconds)\nBear in mind times greater than 20s will take a long time to calculate and the resulting gif will be very large\n");
    scanf("%lg",&temp);
    FILE *thetadata; //Opening connection to output to datafile here, as it will be used by both case options.
    thetadata = fopen("thetadata.txt","w");
    //y mode is the FULL program, including the gnuplot pipe.
    switch(runnable){
        case('y'):
            scalingfactor = 2*(dblpend->l1 + dblpend->l2); //Scaling factor that sets a dynamic range depending on lengths of pendulum chosen

            FILE *grphing; //File for plotting the pendulum

            //Opening the pipe to gnuplot - Setting the output up
            FILE *gnugif = popen("gnuplot -persist","w");
            fprintf(gnugif, "set xrange [%lg:%lg]; set yrange [%lg:%lg]\n",-(scalingfactor),scalingfactor,-(scalingfactor),scalingfactor);
            fprintf(gnugif,"set terminal gif size 600,600 animate delay 1\n");
            fprintf(gnugif,"set output 'pendulumed.gif'\n");
            fprintf(gnugif,"set pointsize 0\n");
            fprintf(gnugif, "set title 'Trajectory of Undamped Double Pendulum'\n");
            fprintf(gnugif, "set style data lines\n");
            //Starting loop - loops for specified amount of time.
            for(i = 0;i<(temp/h);i++){
                //Opening file to contain the coordinates of the simulation - program exits if file cannot be opened.
                if ((grphing=fopen("Coordinates.txt", "w"))==NULL) {
                    fprintf(stdout, "Cannot open Coordinates.txt\n");
                    exit (EXIT_FAILURE);
                }
                //Call to the RK4 method
                rk4(dblpend,h);
                //Setting the coordinates based on rk4 results
                x1 = dblpend->l1*sin(dblpend->t1);
                y1 = -dblpend->l1*cos(dblpend->t1);
                x2 = x1 + dblpend->l2*sin(dblpend->t2);
                y2 = y1 - dblpend->l2*cos(dblpend->t2);
                //prints coordinates to screen for the cool 'hacker' effect.
                printf("%lg\t%lg\t%lg\t%lg\n", x1, y1, x2, y2);
                fprintf(grphing, "%lg\t%lg\t%lg\t%lg\n", x1, y1, x2, y2);
                fprintf(thetadata,"%lg %lg %lg\n",asin(x1/dblpend->l1),asin((x2-x1)/dblpend->l2),t);
                fclose(grphing);



                if (gnugif) {
                    fprintf(gnugif, "unset key; unset ytics; unset xtics\n");
                    fprintf(gnugif, "set multiplot\n");
                    fprintf(gnugif, "plot '-' with lines lw 2 lc rgb 'black', 'Coordinates.txt' u 1:2, 'Coordinates.txt' u 3:4\n");//1:2 being x1,y1.  3:4 is x2,y2.
                    fprintf(gnugif, "%f %f\n", x0, y0);
                    fprintf(gnugif, "%f %f\n", x1, y1);
                    fprintf(gnugif, "%f %f\n", x2, y2);
                    fprintf(gnugif, "e\n");
                    fprintf(gnugif, "\n");
                    fprintf(gnugif, "set nomultiplot\n");
                    fflush(gnugif);
                }
                t += h;

            }
            fclose(thetadata);
            fprintf(gnugif,"exit \n");
            pclose(gnugif);
            printf("You will have a large, unoptimised .gif of the animation of the pendulum in your working directory\nAlso the raw data is there too\n");
            break;

        case('n'):

            //Running the loop but just for data
            for(i = 0;i<(temp/h);i++){
                rk4(dblpend,h);
                x1 = dblpend->l1*sin(dblpend->t1);
                y1 = -dblpend->l1*cos(dblpend->t1);
                x2 = x1 + dblpend->l2*sin(dblpend->t2);
                y2 = y1 - dblpend->l2*cos(dblpend->t2);
                fprintf(thetadata,"%lg %lg %lg\n",asin(x1/dblpend->l1),asin((x2-x1)/dblpend->l2),t);
                t+=h;
            }
            fclose(thetadata);
            printf("Your raw data is stored in t1 t2 t format in the working directory");
            break;
        default:
            printf("Please enter correct input.");
            return 0;
            break;
    }
    printf("\n\nThank you for using the program!!");



    return 0;
}
