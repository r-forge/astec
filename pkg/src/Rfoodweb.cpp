#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <Rcpp.h>
using namespace Rcpp;

void make_average_logbodymass(double *animal_logbodymass, double *animal_bodymass_average, int n_mass){
 if (n_mass>1){
  for (int i=0;i<n_mass;i++){
   animal_bodymass_average[i]=exp(animal_logbodymass[i])*(exp(animal_logbodymass[n_mass])-1)/animal_logbodymass[n_mass];
  }
 }
}

void make_average_logbodymass_matrix(double **animal_logbodymass_matrix, double **animal_bodymass_average_matrix, int *n_mass_list, int n_trophic_group){
 for (int i=0;i<n_trophic_group;i++){
    make_average_logbodymass(animal_logbodymass_matrix[i],animal_bodymass_average_matrix[i],n_mass_list[i]);
 }
}

double animal_respiration(double animal_bodymass, double *metabolic_parameters){
 return metabolic_parameters[0]*exp(metabolic_parameters[1]*log(animal_bodymass));
}

void animal_growth_logbodymass_without_waste(double *animal_density, double *animal_newdensity, double *animal_logbodymass, double *animal_bodymass_average, double *food_assimilated, int n_mass, double *metabolic_parameters){
 for (int i=0;i<n_mass;i++){
  animal_newdensity[i]=animal_density[i];
 }
 for (int i=0;i<n_mass;i++){
  if (animal_density[i]>0.0){
   double respiration=animal_density[i]*animal_respiration((animal_bodymass_average[i]),metabolic_parameters);
   double growing_factor = (food_assimilated[i] - respiration )/ (animal_density[i]*(animal_bodymass_average[i]));
   //waste[i]=((1-animal_growth_parameters[0])*food_assimilated[i])+metabolic_parameters[2]*respiration; // non_assimilated food + part of the respiration that turns into faeces/urine (as opposed to CO2)
   if (growing_factor>=0.0){
    if (i<(n_mass-1)){ // part of the increment goes into the bodymass class above
     double growing_proportion = log(1+growing_factor) / animal_logbodymass[n_mass]; // last cell of bodymass stores the mass width of each cell
     animal_newdensity[i]*=(1.0-growing_proportion)*(1.0+growing_factor);
     animal_newdensity[(i+1)]+= (growing_proportion*(1.0+growing_factor)*animal_density[i]*animal_bodymass_average[i]/animal_bodymass_average[(i+1)]);
    }
    else{ // all the increment stays in the upper bodymass class and is translated as an increase in density
     animal_newdensity[i]*=(1.0+growing_factor);
    }
   }
   else{
    if (i>0){ // part of the decrement goes into the bodymass class below
     double decreasing_proportion = -log(1+growing_factor) / animal_logbodymass[n_mass];
     animal_newdensity[i]*=(1.0-decreasing_proportion)*(1.0+growing_factor);
     animal_newdensity[(i-1)]+= (decreasing_proportion*(1.0+growing_factor)*animal_density[i]*animal_bodymass_average[i]/animal_bodymass_average[(i-1)]);
    }
    else{ // all the decrement stays in the lower bodymass class and is translated as a decrease in density
     animal_newdensity[i]*=(1.0+growing_factor);
    }
   }
  }
 }
}

void animal_growth_logbodymass_matrix_without_waste(double **animal_density_matrix, double **animal_newdensity_matrix, double **animal_logbodymass_matrix, double **animal_bodymass_average_matrix, double **food_assimilated_matrix, int *n_mass_list, double **metabolic_parameters_matrix, int n_trophic_group){
 for (int i=0;i<n_trophic_group;i++){
  animal_growth_logbodymass_without_waste(animal_density_matrix[i],animal_newdensity_matrix[i],animal_logbodymass_matrix[i],animal_bodymass_average_matrix[i],food_assimilated_matrix[i],n_mass_list[i],metabolic_parameters_matrix[i]);
 }
}

void animal_natural_mortality(double *animal_density, double *bodymass_average, double *carrion, int n_mass, double *mortality_parameters, double deltaT){
 for (int i=0;i<n_mass;i++){
    double factor=exp(-mortality_parameters[0]*deltaT);
    carrion[i]=(1.0-factor)*animal_density[i]*bodymass_average[i];
    animal_density[i]*=factor;
 }
}

void animal_natural_mortality_matrix(double **animal_density_matrix, double **bodymass_average_matrix, double **carrion_matrix, int *n_mass_list, double **mortality_parameters_matrix, double deltaT, int n_trophic_group){
 for (int i=0;i<n_trophic_group;i++){
  animal_natural_mortality(animal_density_matrix[i],bodymass_average_matrix[i],carrion_matrix[i],n_mass_list[i],mortality_parameters_matrix[i],deltaT);
 }
}

void animal_reproduction(double *animal_density, double *animal_newdensity, double *animal_bodymass_average, double *juvenile_mass, int *translation_juvenile_mass, double *reproduction_parameters, int n_mass, double deltaT){
 for (int i=0;i<n_mass;i++){
  animal_newdensity[i]=0.0;
 }
 for (int i=0;i<n_mass;i++){
    double juv_mass=animal_density[i]*reproduction_parameters[0]*exp(reproduction_parameters[1]*log(animal_bodymass_average[i]))*deltaT*juvenile_mass[i];
    animal_newdensity[(translation_juvenile_mass[i])]+=(juv_mass/animal_bodymass_average[(translation_juvenile_mass[i])]);
    animal_density[i]-=(juv_mass/animal_bodymass_average[i]);
 }
 for (int i=0;i<n_mass;i++){
  animal_newdensity[i]+=animal_density[i];
 }
}

void animal_reproduction_matrix(double **animal_density_matrix, double **animal_newdensity_matrix, double **animal_bodymass_average_matrix, double **juvenile_mass_matrix, int **translation_juvenile_mass_matrix, double **reproduction_parameters_matrix, int *n_mass_list, double deltaT, int n_trophic_group){
 for (int i=0;i<n_trophic_group;i++){
  animal_reproduction(animal_density_matrix[i],animal_newdensity_matrix[i],animal_bodymass_average_matrix[i],juvenile_mass_matrix[i],translation_juvenile_mass_matrix[i],reproduction_parameters_matrix[i],n_mass_list[i],deltaT);
 }
}

void compute_attack_rate(double **attack_rate_matrix, double **animal_bodymass_average_matrix, double **predation_parameters, int n_trophic_group, int n_carnivore, int *n_mass_list, int *pos_carnivore, int *temp){
 temp[0]=0;
 for (int c=0;c<n_carnivore;c++){
  int i=pos_carnivore[c];
  for (int j=0;j<n_mass_list[i];j++){
   temp[1]=0;
   for (int k=0;k<n_trophic_group;k++){
    for (int l=0;l<n_mass_list[k];l++){
     attack_rate_matrix[(temp[0])][(temp[1])]=predation_parameters[c][0]*animal_bodymass_average_matrix[i][j]*std::exp(-pow((std::log(animal_bodymass_average_matrix[k][l]/animal_bodymass_average_matrix[i][j])-predation_parameters[c][1])/predation_parameters[c][2],2));
     temp[1]+=1;
    }
   }
   temp[0]+=1;
  }
 }
}

void compute_handling_time(double **handling_time_matrix, double **animal_bodymass_average_matrix, double **predation_parameters, int n_trophic_group, int n_carnivore, int *n_mass_list, int *pos_carnivore, int *temp){
    // from Harfoot et al. 2014 - eq.40
 temp[0]=0;
 for (int c=0;c<n_carnivore;c++){
  int i=pos_carnivore[c];
  for (int j=0;j<n_mass_list[i];j++){
   temp[1]=0;
   for (int k=0;k<n_trophic_group;k++){
    for (int l=0;l<n_mass_list[k];l++){
     handling_time_matrix[(temp[0])][(temp[1])]=predation_parameters[c][3]*pow(predation_parameters[c][4]/animal_bodymass_average_matrix[i][j],predation_parameters[c][5])*animal_bodymass_average_matrix[k][l];
     temp[1]+=1;
    }
   }
   temp[0]+=1;
  }
 }
}

void computation_predation_matrix(double **predation_matrix, double *sum_predation_array, int cols_predation_matrix, double **animal_density_matrix,  double **attack_rate_matrix, double **handling_time_matrix, int n_carnivore, int *n_mass_list, int **pos_predation_matrix, int **pos_density_matrix, int *pos_carnivore, int *p, double *temp){
 //computation of F_i,j(t) (eq. 47)
 for (int m=0;m<cols_predation_matrix;m++){
  sum_predation_array[m]=0.0;
 }
 for (int c=0;c<n_carnivore;c++){
  int k=pos_carnivore[c];
  for (int l=0;l<n_mass_list[k];l++){
   p[0]=pos_predation_matrix[k][l];
   //computation of the denominator and storage in temp : 1 + sum_m a_im * d_m * T_im
   temp[0]=1.0;
   for (int m=0;m<cols_predation_matrix;m++){
    temp[0]+=(animal_density_matrix[(pos_density_matrix[m][0])][(pos_density_matrix[m][1])]*attack_rate_matrix[(p[0])][m]*handling_time_matrix[(p[0])][m]);
   }
   //computation of F_i,j(t) and sum_predation_matrix [j] = sum_k F_kj
   for (int m=0;m<cols_predation_matrix;m++){
     predation_matrix[(p[0])][m]=animal_density_matrix[k][l]*attack_rate_matrix[(p[0])][m]/temp[0];
     sum_predation_array[m]+=predation_matrix[(p[0])][m];
   }
  }
 }
}

void animal_predation_with_waste(double **animal_density_matrix, double **animal_consumed_density_matrix, double **animal_bodymass_average_matrix, double **food_assimilated_array, double **waste_matrix, double **predation_matrix, double *sum_predation_array, double deltaT, int n_trophic_species, int n_carnivore, int *n_mass_list, double **animal_growth_parameters_matrix, int **pos_predation_matrix, int *pos_carnivore, int *p, double *temp){
 for (int k=0;k<n_trophic_species;k++){
  for (int l=0;l<n_mass_list[k];l++){
   p[0]=pos_predation_matrix[k][l];
   animal_consumed_density_matrix[k][l]=animal_density_matrix[k][l]*(1-exp(-deltaT*sum_predation_array[(p[0])]));
   animal_density_matrix[k][l]-=animal_consumed_density_matrix[k][l];
   temp[0]=animal_consumed_density_matrix[k][l]*animal_bodymass_average_matrix[k][l]/sum_predation_array[(p[0])]; // biomass amount since multiplication by average bodymass of prey
   // computation of food eaten by each predator (eq. 46)
   for (int c=0;c<n_carnivore;c++){
    int i=pos_carnivore[c];
    for (int j=0;j<n_mass_list[i];j++){
     int m=pos_predation_matrix[i][j];
     food_assimilated_array[m][(p[0])]=temp[0]*predation_matrix[m][(p[0])]; // NB: computation completed two lines below
     waste_matrix[i][j]+=food_assimilated_array[m][(p[0])]*(1-animal_growth_parameters_matrix[i][0]); // NB : wastes need to be set to 0 at each time step.
     food_assimilated_array[m][(p[0])]*=animal_growth_parameters_matrix[i][0];
    }
   }
  }
 }
}

void compile_predator_food_eaten(double **food_assimilated_array,double **food_assimilated_matrix,int lines_predation_matrix,int cols_predation_matrix,int **pos_density_matrix){
 for (int i=0;i<lines_predation_matrix;i++){
  int k=pos_density_matrix[i][0];
  int l=pos_density_matrix[i][1];
  food_assimilated_matrix[k][l]=0.0;
  for (int j=0;j<cols_predation_matrix;j++){
   food_assimilated_matrix[k][l]+=food_assimilated_array[i][j];
  }
 }
}

void nullify(double **waste_matrix,int n_trophic_species,int *n_mass_list){
 for (int i=0;i<n_trophic_species;i++){
    for (int j=0;j<n_mass_list[i];j++){
        waste_matrix[i][j]=0;
    }
 }
}

void collect_detritus(double *detritus,double **waste_matrix,int **translation_waste, double **carrion_matrix, int **translation_carrion,int n_trophic_species,int *n_mass_list){
 for (int i=0;i<n_trophic_species;i++){
    for (int j=0;j<n_mass_list[i];j++){
        detritus[(translation_waste[i][j])]+=waste_matrix[i][j];
        detritus[(translation_carrion[i][j])]+=carrion_matrix[i][j];
    }
 }
}

int spot(double focal_mass, double *array_mass, int n_mass){
 int i=0;
 int res=0;
 while(i<n_mass){
   if (focal_mass>=array_mass[i]){
    res=i;
    i+=1;
   }
   else{
    i=n_mass;
   }
 }
 return res;
}

void compute_herbivory_rate(double **herbivory_rate_matrix,double **animal_bodymass_average_matrix, int *pos_herbivore, double **herbivory_parameters,int n_herbivore, int *n_mass_list, int *p){
 for (int i=0;i<n_herbivore;i++){
  p[0]=pos_herbivore[i];
  for (int j=0;j<n_mass_list[(p[0])];j++){
   herbivory_rate_matrix[i][j]=herbivory_parameters[i][0]*exp(herbivory_parameters[i][1]*log(animal_bodymass_average_matrix[(p[0])][j]));
  }
 }
}

void compute_herbivory_K(double **herbivory_K_matrix,double **animal_bodymass_average_matrix, int *pos_herbivore, double **herbivory_parameters,int n_herbivore, int *n_mass_list, int *p){
 for (int i=0;i<n_herbivore;i++){
  p[0]=pos_herbivore[i];
  for (int j=0;j<n_mass_list[(p[0])];j++){
   herbivory_K_matrix[i][j]=herbivory_parameters[i][2]*animal_bodymass_average_matrix[(p[0])][j];
  }
 }
}

void herbivory_with_waste(double *plant_biomass, double **animal_density_matrix, int *pos_herbivore, double **herbivory_rate_matrix, double **herbivory_K_matrix, double **food_assimilated_matrix, double **animal_growth_parameters, double deltaT, int n_herbivore, int *n_mass_list, double **waste_matrix, double *temp, double *consumed_biomass, int *p){
 temp[0]=0.0;
 for (int i=0;i<n_herbivore;i++){ // computation of F_k
  p[0]=pos_herbivore[i];
  for (int j=0;j<n_mass_list[(p[0])];j++){
   food_assimilated_matrix[(p[0])][j]=herbivory_rate_matrix[i][j]*animal_density_matrix[(p[0])][j]/(plant_biomass[0]+(herbivory_K_matrix[i][j]*animal_density_matrix[(p[0])][j]));
   temp[0]+=food_assimilated_matrix[(p[0])][j]; // computation of sum_k F_k
  }
 }
 consumed_biomass[0]=plant_biomass[0]*(1-exp(-deltaT*temp[0]));
 consumed_biomass[0]/=temp[0]; // for computing speed in herbivory_matrix below
 for (int i=0;i<n_herbivore;i++){ // computation of DeltaB_i amount of plant biomass consumed by each herbivore
  p[0]=pos_herbivore[i];
  for (int j=0;j<n_mass_list[(p[0])];j++){
   food_assimilated_matrix[(p[0])][j]*=consumed_biomass[0];
   waste_matrix[(p[0])][j]=food_assimilated_matrix[(p[0])][j]*(1-animal_growth_parameters[(p[0])][0]);
   food_assimilated_matrix[(p[0])][j]*=animal_growth_parameters[(p[0])][0]; // return the matrix of food assimilated by herbivore
  }
 }
}

void compute_detritivory_rate_matrix3d(double ***detritivory_rate_matrix3d,double **detritivory_rate_matrix, double **detritivory_parameters, double **animal_bodymass_average_matrix, double *detritus, double *detritus_mass, int n_detritus, int n_detritivore, int *pos_detritivore, int *n_mass_list, double *temp, int *p){
 for (int i=0;i<n_detritivore;i++){
  p[0]=pos_detritivore[i];
  for (int j=0;j<n_mass_list[(p[0])];j++){
   temp[0]=0.0;
   for (int k=0;k<n_detritus;k++){
    detritivory_rate_matrix3d[i][j][k]=detritus[k]*detritivory_parameters[i][3]*animal_bodymass_average_matrix[(p[0])][j]*exp(-pow((log(detritus_mass[k]/animal_bodymass_average_matrix[(p[0])][j])-detritivory_parameters[i][4])/detritivory_parameters[i][5],2)); // NB : detritivory has 6 parameters
    temp[0]+=detritivory_rate_matrix3d[i][j][k];
   }
   for (int k=0;k<n_detritus;k++){
    detritivory_rate_matrix3d[i][j][k]*=(detritivory_rate_matrix[i][j]/temp[0]);
   }
  }
 }
}

void detritivory_with_waste(double *detritus, double **animal_density_matrix,int *pos_detritivore, double ***detritivory_rate_matrix3d, double **detritivory_K_matrix, double ***food_assimilated_matrix3d, double **animal_growth_parameters, double deltaT, int n_detritivore, int n_detritus, int *n_mass_list, double **waste_matrix, double *temp, double *consumed_biomass, int *p){
 for (int k=0;k<n_detritus;k++){
  temp[0]=0;
  for (int i=0;i<n_detritivore;i++){ // computation of F_ik
   p[0]=pos_detritivore[i];
   for (int j=0;j<n_mass_list[(p[0])];j++){
    food_assimilated_matrix3d[i][j][k]=detritivory_rate_matrix3d[i][j][k]*animal_density_matrix[(p[0])][j]/(detritus[k]+(detritivory_K_matrix[i][j]*animal_density_matrix[(p[0])][j]));
    temp[0]+=food_assimilated_matrix3d[i][j][k]; // computation of sum_i F_ik
   }
  }
  consumed_biomass[0]=detritus[k]*(1-exp(-deltaT*temp[0]));
  consumed_biomass[0]/=temp[0]; // for computing speed in detritivory_matrix below
  for (int i=0;i<n_detritivore;i++){ // computation of DeltaB_i amount of detritus consumed by each detritivore
   p[0]=pos_detritivore[i];
   for (int j=0;j<n_mass_list[(p[0])];j++){
    food_assimilated_matrix3d[i][j][k]*=consumed_biomass[0];
    waste_matrix[(p[0])][j]+=food_assimilated_matrix3d[i][j][k]*(1-animal_growth_parameters[(p[0])][0]); // waste_matrix needs to be set to zero at each time step
    food_assimilated_matrix3d[i][j][k]*=animal_growth_parameters[(p[0])][0]; // return the matrix of food assimilated by detritivore
   }
  }
 }
}

void compile_detritivore_food_eaten(double **food_assimilated_matrix, double ***food_assimilated_matrix3d, int n_detritivore, int n_detritus, int *pos_detritivore, int *n_mass_list, int *p){
 for (int i=0;i<n_detritivore;i++){
  p[0]=pos_detritivore[i];
  for (int j=0;j<n_mass_list[(p[0])];j++){
   food_assimilated_matrix[(p[0])][j]=0.0;
   for (int k=0;k<n_detritus;k++){
    food_assimilated_matrix[(p[0])][j]+=food_assimilated_matrix3d[i][j][k];
   }
  }
 }
}


bool check(int *trophic_type,int n_carnivore,int n_herbivore,int n_detritivore,int n_trophic_group){
 int t1=0;
 int t2=0;
 int t3=0;
 bool res=false;
 for (int i=0;i<n_trophic_group;i++){
  switch(trophic_type[i]){
   case 0 : t1++;
            break;
   case 1 : t2++;
            break;
   default : t3++;
  }
 }
 if (t1!=n_carnivore){
  res=true;
 }
 if (t2!=n_herbivore){
  res=true;
 }
 if (t3!=n_detritivore){
  res=true;
 }
 if (n_trophic_group!=(n_carnivore+n_herbivore+n_detritivore)){
  res=true;
 }
 for (int i=0;i<n_carnivore;i++){ // carnivores should be in first
  if (trophic_type[i]!=0){
   res=true;
  }
 }
 return res;
}

// [[Rcpp::export]]
List foodweb(List foodweb_inputs){
 //Rcerr <<"astec start"<<"\n";
 /* IMPORTING INPUT FILES */

 int n_trophic_group;
 double deltaT;
 NumericVector glob_param=foodweb_inputs["global_parameters"];

 bool checktest= (glob_param.size()!=6);
 if (checktest){
  List foodweb_outputs = List::create(Named("error") = "invalid global parameters input");
  return foodweb_outputs;
 }


 n_trophic_group=int(glob_param[0]);
 deltaT=glob_param[1];
 deltaT/=365.0; // change of unit: day -> year
 int n_step_dynamics, n_carnivore,n_herbivore,n_detritivore;
 n_step_dynamics=int(glob_param[2]);
 n_carnivore=int(glob_param[3]);
 n_herbivore=int(glob_param[4]);
 n_detritivore=int(glob_param[5]);

 //Rcerr <<"block 1"<<"\n";

 char **name_trophic_groups;
 name_trophic_groups=new char*[(n_trophic_group)];
 List ntg=foodweb_inputs["name_trophic_groups"];

 checktest= (ntg.size()!=n_trophic_group);
 if (checktest){
  List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of global parameters and trophic group names");
  return foodweb_outputs;
 }

 for (int i=0;i<n_trophic_group;i++){
    name_trophic_groups[i]=new char[128];
    String s=ntg[i];
    strcpy(name_trophic_groups[i],s.get_cstring());
 }

 //Rcerr <<"block 2"<<"\n";

 int n_param=26;
 double **trophic_group_parameters;
 trophic_group_parameters= new double*[n_param];
 NumericMatrix tgp=foodweb_inputs["trophic_group_parameters"];
 checktest= ((tgp.nrow()!=n_param)||(tgp.ncol()!=n_trophic_group));
 if (checktest){
  List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of global parameters and trophic group parameters or invalid number of trophic group parameters");
  return foodweb_outputs;
 }

 for (int i=0;i<n_param;i++){
  trophic_group_parameters[i]=new double[n_trophic_group];
  for (int j=0;j<n_trophic_group;j++){
   trophic_group_parameters[i][j]=tgp(i,j);
  }
 }

 //Rcerr <<"block 3"<<"\n";

 int *trophic_type;
 trophic_type=new int[n_trophic_group];
 for (int i=0;i<n_trophic_group;i++){
    trophic_type[i]=int(trophic_group_parameters[3][i]);
 }

 //Rcerr <<"block 4"<<"\n";

 checktest=check(trophic_type,n_carnivore,n_herbivore,n_detritivore,n_trophic_group);
 if (checktest){
  List foodweb_outputs = List::create(Named("error") = "invalid trophic types or problem of ordering : carnivores should be in first");
  return foodweb_outputs;
 }

 //Rcerr <<"block 5"<<"\n";

 int *pos_carnivore;
 pos_carnivore=new int[n_carnivore];
 int *pos_herbivore;
 pos_herbivore=new int[n_herbivore];
 int *pos_detritivore;
 pos_detritivore=new int[n_detritivore];
 int t1=0;
 int t2=0;
 int t3=0;
 for (int i=0;i<n_trophic_group;i++){
    switch(trophic_type[i]){
     case 0 : pos_carnivore[t1]=i;
              t1++;
              break;
     case 1 : pos_herbivore[t2]=i;
              t2++;
              break;
     default : pos_detritivore[t3]=i;
               t3++;
    }
 }

 //Rcerr <<pos_carnivore[0]<<" "<<pos_herbivore[0]<<" "<<pos_detritivore[0]<<" block 6"<<"\n";

 int *n_mass_list;
 n_mass_list= new int[n_trophic_group];
 int cols_predation_matrix=0;
 int lines_predation_matrix=0;
 double **animal_logbodymass_matrix;
 animal_logbodymass_matrix= new double*[n_trophic_group];
 double **animal_bodymass_average_matrix;
 animal_bodymass_average_matrix= new double*[n_trophic_group];
 double **waste_matrix;
 waste_matrix= new double*[n_trophic_group];
 double **carrion_matrix;
 carrion_matrix= new double*[n_trophic_group];
 for (int i=0;i<n_trophic_group;i++){
  n_mass_list[i]=trophic_group_parameters[0][i];
  cols_predation_matrix+=n_mass_list[i];
  if (trophic_type[i]==0){
    lines_predation_matrix+=n_mass_list[i];
  }
  animal_logbodymass_matrix[i]=new double[(n_mass_list[i]+1)];
  if (n_mass_list[i]>1){
   animal_logbodymass_matrix[i][(n_mass_list[i])]=(trophic_group_parameters[2][i]-trophic_group_parameters[1][i])/(n_mass_list[i]-1.0);
  }
  else{
   animal_logbodymass_matrix[i][(n_mass_list[i])]=1.0;
  }
  animal_bodymass_average_matrix[i]=new double[(n_mass_list[i])];
  waste_matrix[i]=new double[(n_mass_list[i])];
  carrion_matrix[i]=new double[(n_mass_list[i])];
  for (int j=0;j<n_mass_list[i];j++){
    animal_logbodymass_matrix[i][j]=trophic_group_parameters[1][i]+j*animal_logbodymass_matrix[i][(n_mass_list[i])];
    animal_bodymass_average_matrix[i][j]=exp(animal_logbodymass_matrix[i][j]); // NB: this will be corrected by the function make_average_logbodymass_matrix below if n_mass_list[i]>1
  }
 }

 //Rcerr <<"block 7"<<"\n";

 double **animal_growth_parameters= new double *[n_trophic_group];
 double **metabolic_parameters= new double *[n_trophic_group];
 double **mortality_parameters= new double *[n_trophic_group];
 double **reproduction_parameters= new double *[n_trophic_group];
 double **predation_parameters;
 predation_parameters= new double *[n_carnivore];
 double **herbivory_parameters;
 herbivory_parameters= new double *[n_herbivore];
 double **detritivory_parameters;
 detritivory_parameters= new double *[n_detritivore];
 t1=0;
 t2=0;
 t3=0;
 for (int i=0;i<n_trophic_group;i++){
  animal_growth_parameters[i]=new double[1];
  metabolic_parameters[i]=new double[3];
  mortality_parameters[i]=new double[1];
  reproduction_parameters[i]=new double[2];
  animal_growth_parameters[i][0]=trophic_group_parameters[4][i];
  for (int j=0;j<3;j++){
   metabolic_parameters[i][j]=trophic_group_parameters[(5+j)][i];
  }
  mortality_parameters[i][0]=trophic_group_parameters[8][i];
  for (int j=0;j<2;j++){
   reproduction_parameters[i][j]=trophic_group_parameters[(9+j)][i];
  }
  switch(trophic_type[i]){
   case 0 : predation_parameters[t1]=new double[6];
            for (int j=0;j<6;j++){
             predation_parameters[t1][j]=trophic_group_parameters[(11+j)][i];
            }
            t1++;
            break;
   case 1 : herbivory_parameters[t2]=new double[3];
            for (int j=0;j<3;j++){
             herbivory_parameters[t2][j]=trophic_group_parameters[(17+j)][i];
            }
            t2++;
            break;
    default : detritivory_parameters[t3]=new double[6];
              for (int j=0;j<6;j++){
               detritivory_parameters[t3][j]=trophic_group_parameters[(20+j)][i];
              }
              t3++;
  }
 }

 Rcerr <<"block 8"<<"\n";

 double **animal_density_matrix;
 animal_density_matrix=new double *[n_trophic_group];
 double **animal_newdensity_matrix;
 animal_newdensity_matrix=new double*[n_trophic_group];
 double **animal_consumed_density_matrix;
 animal_consumed_density_matrix=new double*[n_trophic_group];
 List adm=foodweb_inputs["animal_density_matrix"];
 checktest= (adm.size()!=n_trophic_group);
 if (checktest){
  List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of global parameters and animal density matrix");
  return foodweb_outputs;
 }

 for (int i=0;i<n_trophic_group;i++){
  animal_density_matrix[i]=new double[(n_mass_list[i])];
  animal_newdensity_matrix[i]=new double[(n_mass_list[i])];
  animal_consumed_density_matrix[i]=new double[(n_mass_list[i])];
  NumericVector adm1=adm[i];
  checktest= (adm1.size()!=n_mass_list[i]);
  if (checktest){
   List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of global parameters and animal density matrix");
   return foodweb_outputs;
  }
  for (int j=0;j<n_mass_list[i];j++){
    animal_density_matrix[i][j]=adm1[j];
    animal_newdensity_matrix[i][j]=0;
    animal_consumed_density_matrix[i][j]=0;
  }
 }

 //Rcerr <<"block 9"<<"\n";

 double **juvenile_mass;
 juvenile_mass= new double*[n_trophic_group];
 List jvm=foodweb_inputs["juvenile_mass"];
 checktest= (jvm.size()!=n_trophic_group);
 if (checktest){
  List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of global parameters and juvenile mass");
  return foodweb_outputs;
 }
 for (int i=0;i<n_trophic_group;i++){
  juvenile_mass[i]= new double[(n_mass_list[i])];
  NumericVector jvm1=jvm[i];
  checktest= (jvm1.size()!=n_mass_list[i]);
  if (checktest){
   List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of global parameters and juvenile mass");
   return foodweb_outputs;
  }
  for (int j=0;j<n_mass_list[i];j++){
    juvenile_mass[i][j]=jvm1[j];
  }
 }

 //Rcerr <<"block 10"<<"\n";

 NumericVector dm=foodweb_inputs["detritus_mass"];
 int n_detritus=dm.size();
 double *detritus_mass;
 detritus_mass=new double[n_detritus];
 for (int i=0;i<n_detritus;i++){
    detritus_mass[i]=dm[i];
 }

 Rcerr <<"block 11"<<"\n";

 NumericVector da=foodweb_inputs["detritus_amount"];
 checktest= (da.size()!=n_detritus);
 if (checktest){
  List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of detritus mass and detritus amount");
  return foodweb_outputs;
 }
 double *detritus;
 detritus=new double[n_detritus];
 for (int i=0;i<n_detritus;i++){
    detritus[i]=da[i];
 }

 //Rcerr <<"block 12"<<"\n";

 double **waste_mass;
 waste_mass= new double*[n_trophic_group];
 List wm=foodweb_inputs["waste_mass"];
 checktest= (wm.size()!=n_trophic_group);
 if (checktest){
  List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of global parameters and waste mass");
  return foodweb_outputs;
 }
 for (int i=0;i<n_trophic_group;i++){
  waste_mass[i]= new double[(n_mass_list[i])];
  NumericVector wm1=wm[i];
  checktest= (wm1.size()!=n_mass_list[i]);
  if (checktest){
   List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of global parameters and waste mass");
   return foodweb_outputs;
  }
  for (int j=0;j<n_mass_list[i];j++){
    waste_mass[i][j]=wm1[j];
  }
 }

 //Rcerr <<"block 13"<<"\n";

 double **carrion_mass;
 carrion_mass= new double*[n_trophic_group];
 List cm=foodweb_inputs["carrion_mass"];
 checktest= (cm.size()!=n_trophic_group);
 if (checktest){
  List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of global parameters and carrion mass");
  return foodweb_outputs;
 }
 for (int i=0;i<n_trophic_group;i++){
  carrion_mass[i]= new double[(n_mass_list[i])];
  NumericVector cm1=cm[i];
  checktest= (cm1.size()!=n_mass_list[i]);
  if (checktest){
   List foodweb_outputs = List::create(Named("error") = "inconsistent inputs of global parameters and carrion mass");
   return foodweb_outputs;
  }
  for (int j=0;j<n_mass_list[i];j++){
    carrion_mass[i][j]=cm1[j];
  }
 }

 //Rcerr <<"block 14"<<"\n";

 double *plant_biomass;
 plant_biomass=new double[1]; // no size-structure for plants

 NumericVector plant_biomass_series=foodweb_inputs["plant_biomass_series"];
 int n_plant_biomass_series=plant_biomass_series.size();

 NumericVector plant_senescence_parameters=foodweb_inputs["plant_senescence_parameters"];
 checktest= (plant_senescence_parameters.size()!=2);
 if (checktest){
  List foodweb_outputs = List::create(Named("error") = "inconsistent input of plant senescence parameters");
  return foodweb_outputs;
 }
 double plant_senescence_rate=plant_senescence_parameters[0];
 double mass_plant_litter=plant_senescence_parameters[1];
 int translation_plant_litter=spot(mass_plant_litter,detritus_mass,n_detritus);

 //Rcerr <<"block 15"<<"\n";

 //DEFINE VARIABLES USED IN THE COMPUTATIONS
 int **translation_juvenile;
 translation_juvenile=new int*[n_trophic_group];
 int **translation_waste;
 translation_waste= new int*[n_trophic_group];
 int **translation_carrion;
 translation_carrion= new int*[n_trophic_group];
 for (int i=0;i<n_trophic_group;i++){
  translation_juvenile[i]=new int[(n_mass_list[i])];
  translation_waste[i]=new int[(n_mass_list[i])];
  translation_carrion[i]=new int[(n_mass_list[i])];
  for (int j=0;j<n_mass_list[i];j++){
    translation_juvenile[i][j]=spot(juvenile_mass[i][j],animal_logbodymass_matrix[i],n_mass_list[i]);
    translation_waste[i][j]=spot(waste_mass[i][j],detritus_mass,n_detritus);
    translation_carrion[i][j]=spot(carrion_mass[i][j],detritus_mass,n_detritus);
  }
 }

 int **pos_density_matrix;
 pos_density_matrix=new int*[cols_predation_matrix];
 for (int i=0;i<cols_predation_matrix;i++){
    pos_density_matrix[i]=new int[2];
 }

 int **pos_predation_matrix;
 pos_predation_matrix=new int*[n_trophic_group];
 int pos=0;
 for (int i=0;i<n_trophic_group;i++){
  pos_predation_matrix[i]=new int[(n_mass_list[i])];
  for (int j=0;j<n_mass_list[i];j++){
   pos_predation_matrix[i][j]=pos+j;
   pos_density_matrix[(pos+j)][0]=i;
   pos_density_matrix[(pos+j)][1]=j;
  }
  pos+=n_mass_list[i];
 }

 double **attack_rate_matrix;
 attack_rate_matrix= new double *[lines_predation_matrix];
 double **handling_time_matrix;
 handling_time_matrix= new double *[lines_predation_matrix];
 double **predation_matrix;
 predation_matrix= new double *[lines_predation_matrix];
 double **food_assimilated_array;
 food_assimilated_array=new double*[lines_predation_matrix];
 double *sum_predation_array;
 sum_predation_array= new double[cols_predation_matrix];
 for (int i=0;i<lines_predation_matrix;i++){
  attack_rate_matrix[i]=new double[cols_predation_matrix];
  handling_time_matrix[i]=new double[cols_predation_matrix];
  predation_matrix[i]=new double[cols_predation_matrix];
  food_assimilated_array[i]=new double[cols_predation_matrix];
 }
 double **food_assimilated_matrix;
 food_assimilated_matrix=new double*[n_trophic_group];
 for (int i=0;i<n_trophic_group;i++){
    food_assimilated_matrix[i]=new double[(n_mass_list[i])];
 }

 int *temp;
 temp= new int[2];
 double *temp2;
 temp2= new double[2];
 double *temp3;
 temp3= new double[2];



 // Computing matrices of attack rates and handling times
 // that are independent of population densities
 // and thus only needs to be computed once


 double **herbivory_rate_matrix;
 herbivory_rate_matrix=new double*[n_herbivore];
 double **herbivory_K_matrix;
 herbivory_K_matrix=new double*[n_herbivore];
 for (int i=0;i<n_herbivore;i++){
    herbivory_rate_matrix[i]=new double[(n_mass_list[(pos_herbivore[i])])];
    herbivory_K_matrix[i]=new double[(n_mass_list[(pos_herbivore[i])])];
 }

 double **detritivory_rate_matrix;
 detritivory_rate_matrix=new double*[n_detritivore];
 double **detritivory_K_matrix;
 detritivory_K_matrix=new double*[n_detritivore];
 for (int i=0;i<n_detritivore;i++){
    detritivory_rate_matrix[i]=new double[(n_mass_list[(pos_detritivore[i])])];
    detritivory_K_matrix[i]=new double[(n_mass_list[(pos_detritivore[i])])];
 }
 double ***food_assimilated_matrix3d; // A COMPLETER
 food_assimilated_matrix3d=new double **[n_detritivore];
 double ***detritivory_rate_matrix3d; // A COMPLETER
 detritivory_rate_matrix3d=new double **[n_detritivore];
 for (int i=0;i<n_detritivore;i++){
  food_assimilated_matrix3d[i]=new double*[(n_mass_list[(pos_detritivore[i])])];
  detritivory_rate_matrix3d[i]=new double*[(n_mass_list[(pos_detritivore[i])])];
  for (int j=0;j<(n_mass_list[(pos_detritivore[i])]);j++){
   food_assimilated_matrix3d[i][j]=new double[n_detritus];
   detritivory_rate_matrix3d[i][j]=new double[n_detritus];
  }
 }


 make_average_logbodymass_matrix(animal_logbodymass_matrix,animal_bodymass_average_matrix,n_mass_list,n_trophic_group);
 compute_attack_rate(attack_rate_matrix,animal_bodymass_average_matrix,predation_parameters,n_trophic_group,n_carnivore,n_mass_list,pos_carnivore,temp);
 compute_handling_time(handling_time_matrix,animal_bodymass_average_matrix,predation_parameters,n_trophic_group,n_carnivore,n_mass_list,pos_carnivore,temp);
 compute_herbivory_rate(herbivory_rate_matrix,animal_bodymass_average_matrix,pos_herbivore,herbivory_parameters,n_herbivore,n_mass_list,temp);
 compute_herbivory_K(herbivory_K_matrix,animal_bodymass_average_matrix,pos_herbivore,herbivory_parameters,n_herbivore,n_mass_list,temp);
 compute_herbivory_rate(detritivory_rate_matrix,animal_bodymass_average_matrix,pos_detritivore,detritivory_parameters,n_detritivore,n_mass_list,temp);
 compute_herbivory_K(detritivory_K_matrix,animal_bodymass_average_matrix,pos_detritivore,detritivory_parameters,n_detritivore,n_mass_list,temp);

 //FOOD WEB DYNAMICS
 //predation matrix F_ij


 for (int istep=0;istep<n_step_dynamics;istep++){
  // put waste and carrion to zero
  nullify(waste_matrix,n_trophic_group,n_mass_list);
  nullify(carrion_matrix,n_trophic_group,n_mass_list);

  // read the plant (green) biomass of the corresponding time step.
  plant_biomass[0]=plant_biomass_series[(istep%n_plant_biomass_series)];

  //herbivory and waste computation
  herbivory_with_waste(plant_biomass,animal_density_matrix,pos_herbivore,herbivory_rate_matrix,herbivory_K_matrix,food_assimilated_matrix,animal_growth_parameters,deltaT,n_herbivore,n_mass_list,waste_matrix,temp3,temp2,temp);

  //detritivory and waste computation
  compute_detritivory_rate_matrix3d(detritivory_rate_matrix3d,detritivory_rate_matrix,detritivory_parameters,animal_bodymass_average_matrix,detritus,detritus_mass,n_detritus,n_detritivore,pos_detritivore,n_mass_list,temp2,temp);
  detritivory_with_waste(detritus,animal_density_matrix,pos_detritivore,detritivory_rate_matrix3d,detritivory_K_matrix,food_assimilated_matrix3d,animal_growth_parameters,deltaT,n_detritivore,n_detritus,n_mass_list,waste_matrix,temp3,temp2,temp);
  compile_detritivore_food_eaten(food_assimilated_matrix,food_assimilated_matrix3d,n_detritivore,n_detritus,pos_detritivore,n_mass_list,temp);

  //computation of predation matrix F_ij
  computation_predation_matrix(predation_matrix,sum_predation_array,cols_predation_matrix,animal_density_matrix,attack_rate_matrix,handling_time_matrix,n_carnivore,n_mass_list,pos_predation_matrix,pos_density_matrix,pos_carnivore,temp,temp2);
  //actual predation and food eaten by each size class
  animal_predation_with_waste(animal_density_matrix,animal_consumed_density_matrix,animal_bodymass_average_matrix,food_assimilated_array,waste_matrix,predation_matrix,sum_predation_array,deltaT,n_trophic_group,n_carnivore,n_mass_list,animal_growth_parameters,pos_predation_matrix,pos_carnivore,temp,temp2);
  compile_predator_food_eaten(food_assimilated_array,food_assimilated_matrix,lines_predation_matrix,cols_predation_matrix,pos_density_matrix);
  //growth following food assimilation
  animal_growth_logbodymass_matrix_without_waste(animal_density_matrix,animal_newdensity_matrix,animal_logbodymass_matrix,animal_bodymass_average_matrix,food_assimilated_matrix,n_mass_list,metabolic_parameters,n_trophic_group);
  //natural mortality
  animal_natural_mortality_matrix(animal_newdensity_matrix,animal_bodymass_average_matrix,carrion_matrix,n_mass_list,mortality_parameters,deltaT,n_trophic_group);
  //reproduction
  animal_reproduction_matrix(animal_newdensity_matrix,animal_density_matrix,animal_bodymass_average_matrix,juvenile_mass,translation_juvenile,reproduction_parameters,n_mass_list,deltaT,n_trophic_group);
  //update of the detritus
  collect_detritus(detritus,waste_matrix,translation_waste,carrion_matrix,translation_carrion,n_trophic_group,n_mass_list);
  detritus[translation_plant_litter]+=plant_biomass[0]*plant_senescence_rate; // add plant litter

 }


 Rcerr <<"block 16"<<"\n";

 // convert to List of NumericVector for output (to authorize vectors of different lengths)
 NumericVector ad(n_mass_list[0]);
 for (int i=0;i<n_mass_list[0];i++){
    ad(i)=animal_density_matrix[0][i];
 }
 List list_animal_density= List::create(ad);
 for (int j=1;j<n_trophic_group;j++){
    NumericVector ad2(n_mass_list[j]);
    for (int i=0;i<n_mass_list[j];i++){
     ad2(i)=animal_density_matrix[j][i];
    }
    list_animal_density.push_back(ad2);
 }

 NumericMatrix matrix_attack_rate(lines_predation_matrix,cols_predation_matrix);
 for (int i=0;i<lines_predation_matrix;i++){
  for (int j=0;j<cols_predation_matrix;j++){
    matrix_attack_rate(i,j)=attack_rate_matrix[i][j];
  }
 }

 List foodweb_outputs = List::create(Named("animal_density") = list_animal_density , _["attack_rate"] = matrix_attack_rate);
 return foodweb_outputs;
}

// FUNCTIONS NOT USED ANYMORE

void make_average_bodymass(double *animal_bodymass, double *animal_bodymass_average, int n_mass){
 for (int i=0;i<n_mass;i++){
  animal_bodymass_average[i]=animal_bodymass[i]+0.5*animal_bodymass[n_mass];
 }
}

void make_average_bodymass_matrix(double **animal_bodymass_matrix, double **animal_bodymass_average_matrix, int *n_mass_list, int n_trophic_group){
 for (int i=0;i<n_trophic_group;i++){
    make_average_bodymass(animal_bodymass_matrix[i],animal_bodymass_average_matrix[i],n_mass_list[i]);
 }
}

void animal_growth(double *animal_density, double *animal_newdensity, double *animal_bodymass, double *animal_bodymass_average, double *food_eaten, double *waste, int n_mass, double *animal_growth_parameters, double *metabolic_parameters){
 for (int i=0;i<n_mass;i++){
  animal_newdensity[i]=animal_density[i];
 }
 for (int i=0;i<n_mass;i++){
  if (animal_density[i]>0.0){
   double respiration=animal_density[i]*animal_respiration((animal_bodymass_average[i]),metabolic_parameters);
   double growing_factor = ((animal_growth_parameters[0]*food_eaten[i]) - respiration )/ (animal_density[i]*(animal_bodymass_average[i]));
   waste[i]=((1-animal_growth_parameters[0])*food_eaten[i])+metabolic_parameters[2]*respiration; // non_assimilated food + part of the respiration that turns into faeces/urine (as opposed to CO2)
   if (growing_factor>=0.0){
    if (i<(n_mass-1)){ // part of the increment goes into the bodymass class above
     double growing_threshold = animal_bodymass[(i+1)] / (1.0+growing_factor); // last cell of bodymass stores the mass width of each cell
     double growing_proportion = (animal_bodymass[(i+1)] - growing_threshold) / animal_bodymass[n_mass];
     animal_newdensity[i]*=(1.0-growing_proportion)*(1.0+growing_factor);
     animal_newdensity[(i+1)]+= (growing_proportion*(1.0+growing_factor)*animal_density[i]*animal_bodymass_average[i]/animal_bodymass_average[(i+1)]);
    }
    else{ // all the increment stays in the upper bodymass class and is translated as an increase in density
     animal_newdensity[i]*=(1.0+growing_factor);
    }
   }
   else{
    if (i>0){ // part of the decrement goes into the bodymass class below
     double decreasing_threshold = (animal_bodymass[i]) / (1.0+growing_factor); // last cell of bodymass stores the mass width of each cell
     double decreasing_proportion = (decreasing_threshold-animal_bodymass[i]) / animal_bodymass[n_mass];
     animal_newdensity[i]*=(1.0-decreasing_proportion)*(1.0+growing_factor);
     animal_newdensity[(i-1)]+= (decreasing_proportion*(1.0+growing_factor)*animal_density[i]*animal_bodymass_average[i]/animal_bodymass_average[(i-1)]);
    }
    else{ // all the decrement stays in the lower bodymass class and is translated as a decrease in density
     animal_newdensity[i]*=(1.0+growing_factor);
    }
   }
  }
 }
}

void animal_growth_matrix(double **animal_density_matrix, double **animal_newdensity_matrix, double **animal_bodymass_matrix, double **animal_bodymass_average_matrix, double **food_eaten_matrix, double **waste_matrix, int *n_mass_list, double **animal_growth_parameters_matrix, double **metabolic_parameters_matrix, int n_trophic_group){
 for (int i=0;i<n_trophic_group;i++){
  animal_growth(animal_density_matrix[i],animal_newdensity_matrix[i],animal_bodymass_matrix[i],animal_bodymass_average_matrix[i],food_eaten_matrix[i],waste_matrix[i],n_mass_list[i],animal_growth_parameters_matrix[i],metabolic_parameters_matrix[i]);
 }
}

void animal_growth_logbodymass(double *animal_density, double *animal_newdensity, double *animal_logbodymass, double *animal_bodymass_average, double *food_eaten, double *waste, int n_mass, double *animal_growth_parameters, double *metabolic_parameters){
 for (int i=0;i<n_mass;i++){
  animal_newdensity[i]=animal_density[i];
 }
 for (int i=0;i<n_mass;i++){
  if (animal_density[i]>0.0){
   double respiration=animal_density[i]*animal_respiration((animal_bodymass_average[i]),metabolic_parameters);
   double growing_factor = ((animal_growth_parameters[0]*food_eaten[i]) - respiration )/ (animal_density[i]*(animal_bodymass_average[i]));
   waste[i]=((1-animal_growth_parameters[0])*food_eaten[i])+metabolic_parameters[2]*respiration; // non_assimilated food + part of the respiration that turns into faeces/urine (as opposed to CO2)
   if (growing_factor>=0.0){
    if (i<(n_mass-1)){ // part of the increment goes into the bodymass class above
     double growing_proportion = log(1+growing_factor) / animal_logbodymass[n_mass]; // last cell of bodymass stores the mass width of each cell
     animal_newdensity[i]*=(1.0-growing_proportion)*(1.0+growing_factor);
     animal_newdensity[(i+1)]+= (growing_proportion*(1.0+growing_factor)*animal_density[i]*animal_bodymass_average[i]/animal_bodymass_average[(i+1)]);
    }
    else{ // all the increment stays in the upper bodymass class and is translated as an increase in density
     animal_newdensity[i]*=(1.0+growing_factor);
    }
   }
   else{
    if (i>0){ // part of the decrement goes into the bodymass class below
     double decreasing_proportion = -log(1+growing_factor) / animal_logbodymass[n_mass];
     animal_newdensity[i]*=(1.0-decreasing_proportion)*(1.0+growing_factor);
     animal_newdensity[(i-1)]+= (decreasing_proportion*(1.0+growing_factor)*animal_density[i]*animal_bodymass_average[i]/animal_bodymass_average[(i-1)]);
    }
    else{ // all the decrement stays in the lower bodymass class and is translated as a decrease in density
     animal_newdensity[i]*=(1.0+growing_factor);
    }
   }
  }
 }
}

void animal_growth_logbodymass_matrix(double **animal_density_matrix, double **animal_newdensity_matrix, double **animal_logbodymass_matrix, double **animal_bodymass_average_matrix, double **food_eaten_matrix, double **waste_matrix, int *n_mass_list, double **animal_growth_parameters_matrix, double **metabolic_parameters_matrix, int n_trophic_group){
 for (int i=0;i<n_trophic_group;i++){
  animal_growth_logbodymass(animal_density_matrix[i],animal_newdensity_matrix[i],animal_logbodymass_matrix[i],animal_bodymass_average_matrix[i],food_eaten_matrix[i],waste_matrix[i],n_mass_list[i],animal_growth_parameters_matrix[i],metabolic_parameters_matrix[i]);
 }
}

void animal_predation(double **animal_density_matrix, double **animal_consumed_density_matrix, double **food_eaten_array, double **predation_matrix, double *sum_predation_matrix, int size_predation_matrix, double deltaT, int n_trophic_species, int *n_mass_list, int **pos_predation_matrix, int *p, double *temp){
 for (int k=0;k<n_trophic_species;k++){
  for (int l=0;l<n_mass_list[k];l++){
   p[0]=pos_predation_matrix[k][l];
   animal_consumed_density_matrix[k][l]=predation_matrix[0][(p[0])];
   for (int m=1;m<size_predation_matrix;m++){ // NB: requires that size_predation_matrix>=2
    animal_consumed_density_matrix[k][l]+=predation_matrix[m][(p[0])];
   }
   temp[0]=animal_density_matrix[k][l]*(1-exp(-deltaT*animal_consumed_density_matrix[k][l]));
   animal_consumed_density_matrix[k][l]=temp[0];
   animal_density_matrix[k][l]-=animal_consumed_density_matrix[k][l];
   // computation of food eaten by each predator (eq. 46)
   for (int m=0;m<size_predation_matrix;m++){
    food_eaten_array[m][(p[0])]=temp[0]*predation_matrix[m][(p[0])]/sum_predation_matrix[(p[0])];
   }
  }
 }
}

void read_global_parameters(std::ifstream &read_param, double *global_parameters){
 std::string Text;
 for (int i=0;i<2;i++){
  std::getline(read_param,Text,','); // name of the parameter
  std::getline(read_param,Text);
  global_parameters[i]=std::stoi(Text);
 }
 std::cerr<<"n_trophic_groud read : "<<global_parameters[0]<<std::endl;
 std::cerr<<"deltaT read : "<<global_parameters[1]<<std::endl;
}

void read_trophic_group(std::ifstream &read_param, int n_trophic_group, char **name_trophic_group, double **parameters, int n_param){
 std::string Text;
 // READ NAMES OF TROPHIC GROUPS
 std::getline(read_param,Text,','); // name_trophic_group
 for (int i=0;i<(n_trophic_group-1);i++){
  std::getline(read_param,Text,',');
  strcpy(name_trophic_group[i], Text.c_str());
  std::cerr<<i<<" "<<name_trophic_group[i]<<std::endl;
 }
 std::getline(read_param,Text);
 strcpy(name_trophic_group[(n_trophic_group-1)], Text.c_str());
 std::cerr<<(n_trophic_group-1)<<" "<<name_trophic_group[(n_trophic_group-1)]<<std::endl;

 // READ PARAMETERS OF TROPHIC GROUPS
 for (int i=0;i<n_param;i++){
  std::getline(read_param,Text,','); // name of parameters
  std::cerr<<i<<" : ";
  for (int j=0;j<(n_trophic_group-1);j++){
   std::getline(read_param,Text,',');
   parameters[i][j]=std::stof(Text);
   std::cerr<<parameters[i][j]<<" ";
  }
  std::getline(read_param,Text);
  parameters[i][(n_trophic_group-1)]=std::stof(Text);
  std::cerr<<parameters[i][(n_trophic_group-1)]<<std::endl;
 }
}
