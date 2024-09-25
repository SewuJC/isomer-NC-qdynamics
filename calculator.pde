/*
* Jérôme Chauvet | 25-09-2024,11:50 | France
*
* The following computes and draws the orbits & statistics for the non-commutative dynamics
* of the bistable tautomeric (resp. co-stable conformational isomerization) scheme A <-> B.
*
* This simulation allows you to filter out global thermic extrinsic noise by treating
* the incertainty due to intrinsic non-commutativity and kinetic stochasticity based on 
* temperature dependency independently.
*
* Real-time testings of parameters are done using keyboard commands
*
* REFs : NMR coalescence experiments
*
*/

/* choose representation/graphic output as follows : 
0: time-dependent profile (orbit) ; 1: transversal distribution (statistics) ; 2: P_n+1=f(P_n) phase space ; 3: probability heatmap top-scanning */
float view_mode = 1;    

/* total number of values computated/recorded per orbit (augment for better statistics at the cost of CPU) */
int N = 5000;            

/* TAbles and variables for non-commutative map in first order reversible scheme
---------------------------------------------------------------------------------*/
float[] Ca;
float[] Cb; 
float[] stat_Ca;
float[] stat_Cb;   

/* Matrix A */
float a11;float a12;
float a21;float a22;

/* Matrix B */
float b11;float b12;
float b21;float b22;

/* Matrix AB */
float ab11;float ab12;
float ab21;float ab22;

/* Matrix BA */
float ba11;float ba12;
float ba21;float ba22;

/* value of the chronon and kinetic constants */
float chrn = 0.025 ; // tP = 5.391.10^(-44) s

float ka = 2.56873;// >0
float kb = 2.56873;// >0 

int buff_w = 800;
int buff_h = 300;
PImage buffer = createImage(buff_w, buff_h, RGB);

boolean scan = false;

void setup() {
  size(800, 350);
  background(0);
  strokeWeight(1);
  noSmooth();
  frameRate(50);
  
/* For the A<->B noncommutative map
------------------------------------*/   
  // kinetic constants for each two half-reaction

  // calculate coefficients of matrices A & B 
  a11 = exp(-ka*chrn) ;  a12 = 0;
  a21 = 1-exp(-ka*chrn); a22 = 1;

  b11 = 1 ; b12 = 1-exp(-kb*chrn);
  b21 = 0 ; b22 = exp(-kb*chrn) ;
  
  // create tables for the orbits
  Ca = new float[N];
  Cb = new float[N];
  
  // initiate tables for statistical series
  stat_Ca = new float[width];
  stat_Cb = new float[width];
  
  /* define initial concentrations for product a & product b */
  Ca[0] = 0.8;
  Cb[0] = 0.2;
  
}

void draw() {
  background(0);
   
  String s = "LEFT/RIGHT : -/+ k                UP/DOWN : +/-|ka - kb|                CTRL : new heatmap                SHIFT : view mode";
  fill(200);
  text(s, 20, 20);  
  
  text("ka = "+ka, 20, height-40);text("chronon = "+chrn, 120, height-40);  
  text("kb = "+kb, 20, height-20);
  
  ComputeStepTransform();

  /* draw graphic output 
  ---------------------*/ 
  if (view_mode==0) { // time-dependent profiles (orbits)
    stroke(255);
    for (int i = 1; i < width; i++) {   
      stroke(255,0,0);
      circle(i-1, height-Ca[i-1]*height,1);
      stroke(0,0,255);
      circle(i-1, height-Cb[i-1]*height,1);       
    }
  }

  if (view_mode==1) { // statistical distributions (probabilistic profile)
    stroke(0,0,255);
    int stretchY = 20;
    for (int i = 1; i < width; i++) {   
      stroke(255,0,0);
      line(i-1, height*(1-stat_Ca[i-1]*stretchY), i, height*(1-stat_Ca[i]*stretchY));
      stroke(0,0,255);
      line(i-1, height*(1-stat_Cb[i-1]*stretchY), i, height*(1-stat_Cb[i]*stretchY));
    }
  }

  if (view_mode==2) { // phase space ("i+1 vs i" map)
    stroke(255);
    for (int i = 1; i < width; i++) {
      stroke(255,0,0);
      circle(Ca[i-1]*width, (1-Ca[i])*height, 1);
      stroke(0,0,255);
      circle(Cb[i-1]*width, (1-Cb[i])*height, 1);
    }
  }
  
  if (view_mode==3) { // scan as a heatmap the density of probability over the value of chronon 
    float memorize_chrn = chrn; chrn = 0;    
    if(scan==true){ /* conditional vertical scan whenever CTRL is pressed */    
      buffer.loadPixels();   
        for(int j = 0; j < buff_h; j++){  
          for (int i = 0; i < buff_w; i++) {
            buffer.pixels[i+j*buff_w] = color(floor(255*(1-exp(-stat_Ca[i]*20))),0,floor(255*(1-exp(-stat_Cb[i]*20))));
            //buffer.pixels[i+j*buff_w] = color(0,floor(255*(1-exp(-stat_Ca[i]*20))),0);
          }
          chrn += .001 ;
          ComputeStepTransform();
        }
      buffer.updatePixels();
    }   
    image(buffer, 0, 0);
    scan = false;
    chrn = memorize_chrn ; // refeed original value for chronon   
  } 
}

/* These below are the available keyboard commands and their effect on the simulation
-------------------------------------------------------------------------------------*/
void keyPressed() {
  if (key == CODED) {
    if (keyCode == RIGHT) {
      ka += 0.1;
      kb += 0.1;
      ka=abs(ka);
      kb=abs(kb);
      
      a11 = exp(-ka*chrn) ;  a12 = 0;
      a21 = 1-exp(-ka*chrn); a22 = 1;
      
      b11 = 1 ; b12 = 1-exp(-kb*chrn);
      b21 = 0 ; b22 = exp(-kb*chrn) ;
      
    } if (keyCode == LEFT) {
      ka -= 0.1;
      kb -= 0.1;  
      ka=abs(ka);
      kb=abs(kb);
      
      a11 = exp(-ka*chrn) ;  a12 = 0;
      a21 = 1-exp(-ka*chrn); a22 = 1;
      
      b11 = 1 ; b12 = 1-exp(-kb*chrn);
      b21 = 0 ; b22 = exp(-kb*chrn) ;
      
    } 
    if (keyCode == UP) {
      kb += 0.1;
      kb=abs(kb);
      
      a11 = exp(-ka*chrn) ;  a12 = 0;
      a21 = 1-exp(-ka*chrn); a22 = 1;
      
      b11 = 1 ; b12 = 1-exp(-kb*chrn);
      b21 = 0 ; b22 = exp(-kb*chrn) ;
      
    } if (keyCode == DOWN) {
      kb -= 0.1;
      kb=abs(kb);
      
      a11 = exp(-ka*chrn) ;  a12 = 0;
      a21 = 1-exp(-ka*chrn); a22 = 1;
      
      b11 = 1 ; b12 = 1-exp(-kb*chrn);
      b21 = 0 ; b22 = exp(-kb*chrn) ;     
    }
    
    if (keyCode == CONTROL) { // generate probability heatmap and capture it in buffer
    scan = true;
    }
       
  } if (keyCode == SHIFT) {
    view_mode = (view_mode+1)%4; // when pressed, any other key changes view mode
  }
}

/* This function below computes the dynamics (orbit+statistics) until Nth step 
------------------------------------------------------------------------------*/
void ComputeStepTransform(){  
  for (int i = 1; i < N ; i++) {
      
    float r = random(-1,1);
    
    /* the following models kinetic stochasticity due to ambiant thermic noise */
    float dev_factor = .001;// modulate here standard deviation around value k at will
    float perturb_ka = abs(ka+randomGaussian()*dev_factor);
    float perturb_kb = abs(kb+randomGaussian()*dev_factor);
  
    a11 = exp(-perturb_ka*chrn) ;  a12 = 0;
    a21 = 1-exp(-perturb_ka*chrn); a22 = 1;
  
    b11 = 1 ; b12 = 1-exp(-perturb_kb*chrn);
    b21 = 0 ; b22 = exp(-perturb_kb*chrn) ;
    
    /* reduce incertainty due to intrinsic non-commutativity by choosing AB or BA at random */
    if(r>0){    
      Ca[i] = (a11+(1-b22)*(1-a11))*Ca[i-1] + (1-b22)*Cb[i-1] ;
      Cb[i] = b22*(1-a11)          *Ca[i-1] + b22    *Cb[i-1];    
    }
    if(r<=0){       
       Ca[i] = a11    *Ca[i-1] + a11*(1-b22)          *Cb[i-1];
       Cb[i] = (1-a11)*Ca[i-1] + ((1-a11)*(1-b22)+b22)*Cb[i-1];      
    }
  }

  /* collect statistical points from the values of orb() */
  for (int i = 1; i < N; i++) {
    stat_Ca[floor(Ca[i]*(width-1))] +=1 ;
    stat_Cb[floor(Cb[i]*(width-1))] +=1 ;
  }
  
  /* normalize statistical values */
  for (int i = 0; i < width; i++) {
    stat_Ca[i] /= N;
    stat_Cb[i] /= N;
  }   
}
