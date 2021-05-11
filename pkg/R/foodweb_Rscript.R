#########################
## internal functions
#########################

Rmake_average_logbodymass<-function(animal_logbodymass,n_mass){
 if (n_mass>1){
  return(exp(animal_logbodymass[1:n_mass])*(exp(animal_logbodymass[(n_mass+1)])-1)/animal_logbodymass[(n_mass+1)])
 }
 else{
  return(exp(animal_logbodymass))
 }
}

Rmake_average_logbodymass_list<-function(animal_logbodymass_list,n_mass_list,n_trophic_group){
 res=vector("list", n_trophic_group)
 for (i in 1:n_trophic_group){
    res[[i]]=Rmake_average_logbodymass(animal_logbodymass_list[[i]],n_mass_list[i])
 }
 res
}

Ranimal_respiration<-function(animal_bodymass_average,metabolic_parameters,deltaT){
 deltaT*metabolic_parameters[1]*exp(metabolic_parameters[2]*log(animal_bodymass_average))
}

Ranimal_growth_logbodymass_without_waste<-function(animal_density,animal_logbodymass,animal_bodymass_average,food_assimilated,n_mass,metabolic_parameters,deltaT){
 res=animal_density
 changing_size=res*0.0
 respiration=animal_density*Ranimal_respiration(animal_bodymass_average,metabolic_parameters,deltaT)
 for (i in 1:n_mass){
  if (animal_density[i]>0.0){
   growing_factor = max(-1.0,((food_assimilated[i] - respiration[i] )/ (animal_density[i]*animal_bodymass_average[i]))) # organisms cannot respire more than their body mass.
   if (growing_factor>=0.0){
    if (i<n_mass){ # part of the increment goes into the bodymass class above
     growing_proportion= min(1.0,(log(1.0+growing_factor) / animal_logbodymass[(n_mass+1)])) # last cell of bodymass stores the mass width of each cell
     res[i]=(res[i]*(1.0-growing_proportion)*(1.0+growing_factor))
     changing_size[(i+1)]=(changing_size[(i+1)]+(growing_proportion*(1.0+growing_factor)*animal_density[i]*animal_bodymass_average[i]/animal_bodymass_average[(i+1)]))
    }
    else{ # all the increment stays in the upper bodymass class and is translated as an increase in density
     res[i]=(res[i]*(1.0+growing_factor))
    }
   }
   else{
    if (i>1){ # part of the decrement goes into the bodymass class below
     decreasing_proportion = min(1.0,(-log(1.0+growing_factor) / animal_logbodymass[(n_mass+1)]))
     res[i]=(res[i]*(1.0-decreasing_proportion)*(1.0+growing_factor))
     changing_size[(i-1)]= (changing_size[(i-1)]+(decreasing_proportion*(1.0+growing_factor)*animal_density[i]*animal_bodymass_average[i]/animal_bodymass_average[(i-1)]))
    }
    else{ # all the decrement stays in the lower bodymass class and is translated as a decrease in density
     res[i]=(res[i]*(1.0+growing_factor))
    }
   }
  }
 }
 res+changing_size
}

Ranimal_growth_logbodymass_matrix_without_waste<-function(animal_density_list,animal_logbodymass_list,animal_bodymass_average_list,food_assimilated_list,n_mass_list,metabolic_parameters_matrix,n_trophic_group,deltaT){
 res=vector("list",n_trophic_group)
 for (i in 1:n_trophic_group){
  res[[i]]=Ranimal_growth_logbodymass_without_waste(animal_density_list[[i]],animal_logbodymass_list[[i]],animal_bodymass_average_list[[i]],food_assimilated_list[[i]],n_mass_list[i],metabolic_parameters_matrix[i,],deltaT)
 }
 res
}

Ranimal_natural_mortality<-function(animal_density,bodymass_average,n_mass,mortality_parameters,deltaT){
 factor=exp(-deltaT*mortality_parameters[1]*exp(mortality_parameters[2]*log(bodymass_average)))
 carrion=(1.0-factor)*animal_density*bodymass_average
 newdens=animal_density*factor
 res=list(newdens,carrion)
 res
}

Ranimal_natural_mortality_matrix<-function(animal_density_list,bodymass_average_list,n_mass_list,mortality_parameters_matrix,deltaT,n_trophic_group){
 res=vector("list",2)
 res[[1]]=vector("list",n_trophic_group)
 res[[2]]=vector("list",n_trophic_group)
 for (i in 1:n_trophic_group){
  temp=Ranimal_natural_mortality(animal_density_list[[i]],bodymass_average_list[[i]],n_mass_list[i],mortality_parameters_matrix[i,],deltaT)
  res[[1]][[i]]=temp[[1]]
  res[[2]][[i]]=temp[[2]]
 }
  res
}

Ranimal_reproduction<-function(animal_density,animal_bodymass_average,juvenile_mass,translation_juvenile_mass,reproduction_parameters,n_mass,deltaT){
 animal_newdens=array(0,n_mass)
 animal_dens_adult=animal_density
 juv_mass=animal_density*reproduction_parameters[1]*exp(reproduction_parameters[2]*log(animal_bodymass_average))*deltaT*juvenile_mass
 for (i in 1:n_mass){
  animal_newdens[(translation_juvenile_mass[i])]=animal_newdens[(translation_juvenile_mass[i])]+(juv_mass[i]/animal_bodymass_average[(translation_juvenile_mass[i])])
  animal_dens_adult[i]=animal_dens_adult[i]-(juv_mass[i]/animal_bodymass_average[i])
 }
 animal_newdens+animal_dens_adult
}

Ranimal_reproduction_matrix<-function(animal_density_list,animal_bodymass_average_list,juvenile_mass_list,translation_juvenile_mass_list,reproduction_parameters_matrix,n_mass_list,deltaT,n_trophic_group){
 res=vector("list",n_trophic_group)
 for (i in 1:n_trophic_group){
  res[[i]]=Ranimal_reproduction(animal_density_list[[i]],animal_bodymass_average_list[[i]],juvenile_mass_list[[i]],translation_juvenile_mass_list[[i]],reproduction_parameters_matrix[i,],n_mass_list[i],deltaT)
 }
 res
}

Rcompute_attack_rate<-function(animal_bodymass_average_list,predation_parameters,n_trophic_group,n_carnivore,n_mass_list,pos_carnivore){
 attack_rate_matrix=matrix(0,sum(n_mass_list[pos_carnivore]),sum(n_mass_list))
 p1=1
 for (c in 1:n_carnivore){
  i=pos_carnivore[c]
  for (j in 1:n_mass_list[i]){
   p2=1
   for (k in 1:n_trophic_group){
    for (l in 1:n_mass_list[k]){
     attack_rate_matrix[p1,p2]=predation_parameters[c,1]*animal_bodymass_average_list[[i]][j]*exp(-((log(animal_bodymass_average_list[[k]][l]/animal_bodymass_average_list[[i]][j])-predation_parameters[c,2])/predation_parameters[c,3])^2)
     p2=p2+1
    }
   }
   p1=p1+1
  }
 }
 attack_rate_matrix
}

Rcompute_handling_time<-function(animal_bodymass_average_list,predation_parameters,n_trophic_group,n_carnivore,n_mass_list,pos_carnivore){
    # from Harfoot et al. 2014 - eq.40
 handling_time_matrix=matrix(0,sum(n_mass_list[pos_carnivore]),sum(n_mass_list))
 p1=1
 for (c in 1:n_carnivore){
  i=pos_carnivore[c]
  for (j in 1:n_mass_list[i]){
   p2=1
   for (k in 1:n_trophic_group){
    for (l in 1:n_mass_list[k]){
     handling_time_matrix[p1,p2]=predation_parameters[c,4]*((predation_parameters[c,5]/animal_bodymass_average_list[[i]][j])^predation_parameters[c,6])*animal_bodymass_average_list[[k]][l]
     p2=p2+1
    }
   }
   p1=p1+1
  }
 }
 handling_time_matrix
}

Rcomputation_predation_matrix<-function(cols_predation_matrix,animal_density_list,attack_rate_matrix,handling_time_matrix,n_carnivore,n_mass_list,pos_density_matrix,pos_carnivore){
 # computation of F_i,j(t) (eq. 47)
 predation_matrix=matrix(0,sum(n_mass_list[pos_carnivore]),sum(n_mass_list))
 sum_predation_array=array(0,cols_predation_matrix)
 p1=1
 for (c in 1:n_carnivore){
  k=pos_carnivore[c]
  for (l in 1:n_mass_list[k]){
   # computation of the denominator and storage in temp : 1 + sum_m a_im * d_m * T_im
   temp=1.0
   for (m in 1:cols_predation_matrix){
    temp=temp+(animal_density_list[[(pos_density_matrix[m,1])]][(pos_density_matrix[m,2])]*attack_rate_matrix[p1,m]*handling_time_matrix[p1,m])
   }
   # computation of F_i,j(t) and sum_predation_matrix [j] = sum_k F_kj
   for (m in 1:cols_predation_matrix){
     predation_matrix[p1,m]=animal_density_list[[k]][l]*attack_rate_matrix[p1,m]/temp
     sum_predation_array[m]=(sum_predation_array[m]+predation_matrix[p1,m])
   }
   p1=p1+1
  }
 }
 list(predation_matrix,sum_predation_array)
}

Ranimal_predation_with_waste<-function(waste_list,animal_density_list,animal_bodymass_average_list,predation_matrix,sum_predation_array,deltaT,n_trophic_species,n_carnivore,n_mass_list,animal_growth_parameters_matrix,pos_predation_list,pos_carnivore){
 newwaste_list=waste_list
 animal_consumed_density_list=vector("list",n_trophic_species)
 animal_newdensity_list=animal_density_list
 food_assimilated_matrix=matrix(0,sum(n_mass_list),sum(n_mass_list))
 for (k in 1:n_trophic_species){
  animal_consumed_density_list[[k]]=array(0,n_mass_list[k])
 }
 for (k in 1:n_trophic_species){
  for (l in 1:n_mass_list[k]){
   p1=pos_predation_list[[k]][l]
   animal_consumed_density_list[[k]][l]=animal_density_list[[k]][l]*(1-exp(-deltaT*sum_predation_array[p1]))
   animal_newdensity_list[[k]][l]=animal_density_list[[k]][l]-animal_consumed_density_list[[k]][l]
   temp=animal_consumed_density_list[[k]][l]*animal_bodymass_average_list[[k]][l]
   if (sum_predation_array[p1]>0.0){
    temp=(temp/sum_predation_array[p1]) # biomass amount since multiplication by average bodymass of prey
   }
   # computation of food eaten by each predator (eq. 46)
   mc=1
   for (c in 1:n_carnivore){
    i=pos_carnivore[c]
    for (j in 1:n_mass_list[i]){
     m=pos_predation_list[[i]][j]
     food_assimilated_matrix[m,p1]=temp*predation_matrix[mc,p1] # NB: computation completed two lines below
     newwaste_list[[i]][j]=newwaste_list[[i]][j]+(food_assimilated_matrix[m,p1]*(1-animal_growth_parameters_matrix[i,1])) # NB : wastes need to be set to 0 at each time step.
     food_assimilated_matrix[m,p1]=(food_assimilated_matrix[m,p1]*animal_growth_parameters_matrix[i,1])
     mc=mc+1
    }
   }
  }
 }
 list(animal_consumed_density_list,animal_newdensity_list,food_assimilated_matrix,newwaste_list)
}

Rcompile_predator_food_eaten<-function(food_assimilated_array,food_assimilated_list,cols_predation_matrix,pos_predation_list,n_carnivore,pos_carnivore,n_mass_list){
 newfood_assimilated_list=food_assimilated_list
 for (c in 1:n_carnivore){
  i=pos_carnivore[c]
  for (j in 1:n_mass_list[i]){
   newfood_assimilated_list[[i]][j]=0.0
   l=pos_predation_list[[i]][j]
   for (k in 1:cols_predation_matrix){
    newfood_assimilated_list[[i]][j]=(newfood_assimilated_list[[i]][j]+food_assimilated_array[l,k])
   }
  }
 }
 newfood_assimilated_list
}

Rnullify<-function(n_trophic_species,n_mass_list){
 new_waste_list=vector("list",n_trophic_species)
 for (i in 1:n_trophic_species){
  for (j in 1:n_mass_list[i]){
   new_waste_list[[i]][j]=0
  }
 }
 new_waste_list
}

Rnullify2<-function(n_trophic_species,n_mass_list,n_det){
 new_waste_list=vector("list",n_trophic_species)
 for (i in 1:n_trophic_species){
  new_waste_list[[i]]=matrix(0,n_mass_list[i],n_det)
 }
 new_waste_list
}

Rcollect_detritus<-function(detritus,waste_list,translation_waste,carrion_list,translation_carrion,n_trophic_species,n_mass_list){
 newdetritus=detritus
 for (i in 1:n_trophic_species){
    for (j in 1:n_mass_list[i]){
        newdetritus[(translation_waste[[i]][j])]=(newdetritus[(translation_waste[[i]][j])]+waste_list[[i]][j])
        newdetritus[(translation_carrion[[i]][j])]=(newdetritus[(translation_carrion[[i]][j])]+carrion_list[[i]][j])
    }
 }
 newdetritus
}

Rspot<-function(focal_mass,array_mass,n_mass){
 i=1
 res=1
 while(i<=n_mass){
   if (focal_mass>=array_mass[i]){
    res=i
    i=i+1
   }
   else{
    i=n_mass+1
   }
 }
 res
}

Rcompute_herbivory_rate<-function(animal_bodymass_average_list,pos_herbivore,herbivory_parameters,n_herbivore,n_mass_list){
 herbivory_rate_list=vector("list",n_herbivore)
 for (i in 1:n_herbivore){
  p1=pos_herbivore[i]
  herbivory_rate_list[[i]]=array(0,n_mass_list[p1])
  for (j in 1:n_mass_list[p1]){
   herbivory_rate_list[[i]][j]=herbivory_parameters[i,1]*exp(herbivory_parameters[i,2]*log(animal_bodymass_average_list[[p1]][j]))
  }
 }
 herbivory_rate_list
}

Rcompute_herbivory_K<-function(animal_bodymass_average_list,pos_herbivore,herbivory_parameters,n_herbivore,n_mass_list){
 herbivory_K_list=vector("list",n_herbivore)
 for (i in 1:n_herbivore){
  p1=pos_herbivore[i]
  herbivory_K_list[[i]]=array(0,n_mass_list[p1])
  for (j in 1:n_mass_list[p1]){
   herbivory_K_list[[i]][j]=herbivory_parameters[i,3]*exp(herbivory_parameters[i,4]*log(animal_bodymass_average_list[[p1]][j]))
  }
 }
 herbivory_K_list
}

Rcompute_detritivory_rates<-function(animal_bodymass_average_list,pos_detritivore,detritivory_parameters,n_detritivore,n_mass_list){
 detritivory_attack_rate=vector("list",n_detritivore)
 detritivory_handling_time=vector("list",n_detritivore)
 for (i in 1:n_detritivore){
  detritivory_attack_rate[[i]]=array(0,n_mass_list[(pos_detritivore[i])])
  detritivory_handling_time[[i]]=array(0,n_mass_list[(pos_detritivore[i])])
 }
 for (i in 1:n_detritivore){
  p1=pos_detritivore[i]
  for (j in 1:n_mass_list[p1]){
   detritivory_attack_rate[[i]][j]=detritivory_parameters[i,1]*exp(detritivory_parameters[i,2]*log(animal_bodymass_average_list[[p1]][j]))
   detritivory_handling_time[[i]][j]=detritivory_parameters[i,3]*exp(detritivory_parameters[i,4]*log(animal_bodymass_average_list[[p1]][j]))
  }
 }
 list(detritivory_attack_rate,detritivory_handling_time)
}

Rherbivory_with_waste<-function(plant_biomass,animal_density_list,pos_herbivore,herbivory_rate_list,herbivory_K_list,food_assimilated_list,animal_growth_parameters,deltaT,n_herbivore,n_mass_list, waste_list){
 newfood_assimilated_list=food_assimilated_list
 newwaste_list=waste_list
 temp=0.0
 for (i in 1:n_herbivore){ # computation of F_k
  p1=pos_herbivore[i]
  for (j in 1:n_mass_list[p1]){
   newfood_assimilated_list[[p1]][j]=herbivory_rate_list[[i]][j]*animal_density_list[[p1]][j]/(plant_biomass[1]+(herbivory_K_list[[i]][j]*animal_density_list[[p1]][j]))
   temp=(temp+newfood_assimilated_list[[p1]][j]) # computation of sum_k F_k
  }
 }
 consumed_biomass=plant_biomass[1]*(1-exp(-deltaT*temp))
 newplant_biomass=plant_biomass-consumed_biomass
 if (temp>0.0){
  consumed_biomass=(consumed_biomass/temp) # for computing speed in herbivory_matrix below
 }
 for (i in 1:n_herbivore){ # computation of DeltaB_i amount of plant biomass consumed by each herbivore
  p1=pos_herbivore[i]
  for (j in 1:n_mass_list[p1]){
   newfood_assimilated_list[[p1]][j]=(newfood_assimilated_list[[p1]][j]*consumed_biomass)
   newwaste_list[[p1]][j]=newfood_assimilated_list[[p1]][j]*(1-animal_growth_parameters[p1,1])
   newfood_assimilated_list[[p1]][j]=(newfood_assimilated_list[[p1]][j]*animal_growth_parameters[p1,1]) # return the matrix of food assimilated by herbivore
  }
 }
 list(newfood_assimilated_list,newwaste_list,newplant_biomass)
}

Rdetritivory_with_waste<-function(detritus,animal_density_list,pos_detritivore,detritivory_attack_rate,detritivory_handling_time,food_assimilated_matrix3d,animal_growth_parameters,deltaT,n_detritivore,n_detritus,n_mass_list,waste_list){
 newfood_assimilated_matrix3d=food_assimilated_matrix3d
 newwaste_list=waste_list
 newdetritus=detritus
 for (k in 1:n_detritus){
  temp=0
  for (i in 1:n_detritivore){ # computation of F_ik
   p1=pos_detritivore[i]
   for (j in 1:n_mass_list[p1]){
    newfood_assimilated_matrix3d[[i]][j,k]=detritivory_attack_rate[[i]][j]*animal_density_list[[p1]][j]*detritus[k]/(1+detritus[k]*detritivory_attack_rate[[i]][j]*detritivory_handling_time[[i]][j])
    temp=(temp+newfood_assimilated_matrix3d[[i]][j,k]) # computation of sum_i F_ik
   }
  }
  consumed_biomass=detritus[k]*(1-exp(-deltaT*temp))
  newdetritus[k]=detritus[k]-consumed_biomass
  if (temp>0.0){
   consumed_biomass=(consumed_biomass/temp) # for computing speed in detritivory_matrix below
  }
  for (i in 1:n_detritivore){ # computation of DeltaB_i amount of detritus consumed by each detritivore
   p1=pos_detritivore[i]
   for (j in 1:n_mass_list[p1]){
    newfood_assimilated_matrix3d[[i]][j,k]=(newfood_assimilated_matrix3d[[i]][j,k]*consumed_biomass)
    newwaste_list[[p1]][j]=(newwaste_list[[p1]][j]+newfood_assimilated_matrix3d[[i]][j,k]*(1-animal_growth_parameters[p1,1])) # waste_matrix needs to be set to zero at each time step
    newfood_assimilated_matrix3d[[i]][j,k]=(newfood_assimilated_matrix3d[[i]][j,k]*animal_growth_parameters[p1,1]) # return the matrix of food assimilated by detritivore
   }
  }
 }
 list(newfood_assimilated_matrix3d,newwaste_list,newdetritus)
}

Rcompile_detritivore_food_eaten<-function(food_assimilated_list,food_assimilated_matrix3d,n_detritivore,n_detritus,pos_detritivore,n_mass_list){
 newfood_assimilated_list=food_assimilated_list
 for (i in 1:n_detritivore){
  p1=pos_detritivore[i]
  for (j in 1:n_mass_list[p1]){
   newfood_assimilated_list[[p1]][j]=0.0
   for (k in 1:n_detritus){
    newfood_assimilated_list[[p1]][j]=(newfood_assimilated_list[[p1]][j]+food_assimilated_matrix3d[[i]][j,k])
   }
  }
 }
 newfood_assimilated_list
}

Rcheck<-function(trophic_type,n_carnivore,n_herbivore,n_detritivore,n_trophic_group){
 t1=0
 t2=0
 t3=0
 res=FALSE
 for (i in 1:n_trophic_group){
  if(trophic_type[i]==0){
   t1=t1+1
  }
  else{
   if (trophic_type[i]==1){
    t2=t2+1
   }
   else{
    if (trophic_type[i]==2){
     t3=t3+1
    }
   }
  }
 }
 if (t1!=n_carnivore){
  res=TRUE
 }
 if (t2!=n_herbivore){
  res=TRUE
 }
 if (t3!=n_detritivore){
  res=TRUE
 }
 if (n_trophic_group!=(n_carnivore+n_herbivore+n_detritivore)){
  res=TRUE
 }
 for (i in 1:n_carnivore){ # carnivores should be in first
  if (trophic_type[i]!=0){
   res=TRUE
  }
 }
 res
}


 #####################
 #####################
 ### MAIN FUNCTION ###
 #####################
 #####################

Rfoodweb<-function(foodweb_inputs){

 ##########################################
 ### IMPORTING AND CHECKING foodweb_inputs
 ##########################################

 glob_param=foodweb_inputs$global_parameters
 checktest= (length(glob_param)!=7)
 if (checktest){
  err="invalid global parameters input"
  return(err)
 }

 n_trophic_group=glob_param[1]
 deltaT=glob_param[2] # deltaT expressed in days
 n_step_dynamics=glob_param[3]
 n_carnivore=glob_param[4]
 n_herbivore=glob_param[5]
 n_detritivore=glob_param[6]
 n_step_output=glob_param[7]

 name_trophic_groups=foodweb_inputs$name_trophic_groups
 checktest= (length(name_trophic_groups)!=n_trophic_group)
 if (checktest){
  err="inconsistent inputs of global parameters and trophic group names"
  return(err)
 }

 n_param=26
 trophic_group_parameters=foodweb_inputs$trophic_group_parameters
 checktest= ((dim(trophic_group_parameters)[1]!=n_param)||(dim(trophic_group_parameters)[2]!=n_trophic_group))
 if (checktest){
  err="inconsistent inputs of global parameters and trophic group parameters or invalid number of trophic group parameters"
  return(err)
 }

 trophic_type=trophic_group_parameters[4,]
 checktest=Rcheck(trophic_type,n_carnivore,n_herbivore,n_detritivore,n_trophic_group)
 if (checktest){
  err="invalid trophic types or problem of ordering : carnivores should be in first"
  return(err)
 }


 pos_carnivore=array(0,n_carnivore)
 pos_herbivore=array(0,n_herbivore)
 pos_detritivore=array(0,n_detritivore)
 t1=1
 t2=1
 t3=1
 for (i in 1:n_trophic_group){
  if (trophic_type[i]==0){
   pos_carnivore[t1]=i
   t1=t1+1
  }
  if (trophic_type[i]==1){
   pos_herbivore[t2]=i
   t2=t2+1
  }
  if (trophic_type[i]==2){
   pos_detritivore[t3]=i
   t3=t3+1
  }
 }

 n_mass_list=array(0,n_trophic_group)
 animal_logbodymass_list=vector("list",n_trophic_group)
 animal_bodymass_average_list=vector("list",n_trophic_group)
 waste_list=vector("list",n_trophic_group)
 carrion_list=vector("list",n_trophic_group)
 for (i in 1:n_trophic_group){
  n_mass_list[i]=trophic_group_parameters[1,i]
  animal_logbodymass_list[[i]]=array(0,(n_mass_list[i]+1))
  if (n_mass_list[i]>1){
   animal_logbodymass_list[[i]][(n_mass_list[i]+1)]=(trophic_group_parameters[3,i]-trophic_group_parameters[2,i])/(n_mass_list[i]-1.0)
  }
  else{
   animal_logbodymass_list[[i]][(n_mass_list[i]+1)]=1.0
  }
  animal_bodymass_average_list[[i]]=array(0,n_mass_list[i])
  waste_list[[i]]=array(0,n_mass_list[i])
  carrion_list[[i]]=array(0,n_mass_list[i])
  for (j in 1:n_mass_list[i]){
    animal_logbodymass_list[[i]][j]=trophic_group_parameters[2,i]+(j-1)*animal_logbodymass_list[[i]][(n_mass_list[i]+1)]
    animal_bodymass_average_list[[i]][j]=exp(animal_logbodymass_list[[i]][j]) # NB: this will be corrected by the function make_average_logbodymass_matrix below if n_mass_list[i]>1
  }
 }

 cols_predation_matrix=sum(n_mass_list)
 lines_predation_matrix=sum(n_mass_list[pos_carnivore])


 animal_bodymass_average_list=Rmake_average_logbodymass_list(animal_logbodymass_list,n_mass_list,n_trophic_group)

 animal_growth_parameters=matrix(0,n_trophic_group,1)
 metabolic_parameters=matrix(0,n_trophic_group,3)
 mortality_parameters=matrix(0,n_trophic_group,2)
 reproduction_parameters=matrix(0,n_trophic_group,2)
 predation_parameters=matrix(0,n_carnivore,6)
 herbivory_parameters=matrix(0,n_herbivore,4)
 detritivory_parameters=matrix(0,n_detritivore,4)
 t1=1
 t2=1
 t3=1
 for (i in 1:n_trophic_group){
  animal_growth_parameters[i,1]=trophic_group_parameters[5,i]
  for (j in 1:3){
   metabolic_parameters[i,j]=trophic_group_parameters[(5+j),i]
  }
  mortality_parameters[i,1]=trophic_group_parameters[9,i]
  mortality_parameters[i,2]=trophic_group_parameters[10,i]
  for (j in 1:2){
   reproduction_parameters[i,j]=trophic_group_parameters[(10+j),i]
  }
  if (trophic_type[i]==0){
   for (j in 1:6){
    predation_parameters[t1,j]=trophic_group_parameters[(12+j),i]
   }
   t1=t1+1
  }
  if (trophic_type[i]==1){
   for (j in 1:4){
    herbivory_parameters[t2,j]=trophic_group_parameters[(18+j),i]
   }
   t2=t2+1
  }
  if (trophic_type[i]==2){
   for (j in 1:4){
    detritivory_parameters[t3,j]=trophic_group_parameters[(22+j),i]
   }
   t3=t3+1
  }
 }

 animal_density_list=vector("list",n_trophic_group)
 adm=foodweb_inputs$animal_density_matrix
 checktest= (length(adm)!=n_trophic_group)
 if (checktest){
  err="inconsistent inputs of global parameters and animal density matrix"
  return(err)
 }
 for (i in 1:n_trophic_group){
  checktest= (length(adm[[i]])!=n_mass_list[i])
  if (checktest){
   err="inconsistent inputs of global parameters and animal density matrix"
   return(err)
  }
  animal_density_list[[i]]=adm[[i]]
 }

 juvenile_mass_list=foodweb_inputs$juvenile_mass
 checktest= (length(juvenile_mass_list)!=n_trophic_group)
 if (checktest){
  err="inconsistent inputs of global parameters and juvenile mass"
  return(err)
 }
 for (i in 1:n_trophic_group){
  checktest= (length(juvenile_mass_list[[i]])!=n_mass_list[i])
  if (checktest){
   err="inconsistent inputs of global parameters and juvenile mass"
   return(err)
  }
 }

 detritus_mass=foodweb_inputs$detritus_mass
 n_detritus=length(detritus_mass)

 detritus=foodweb_inputs$detritus_amount
 checktest= (length(detritus)!=n_detritus)
 if (checktest){
  err="inconsistent inputs of detritus mass and detritus amount"
  return(err)
 }

 waste_mass=foodweb_inputs$waste_mass
 checktest= (length(waste_mass)!=n_trophic_group)
 if (checktest){
  err="inconsistent inputs of detritus mass and detritus amount"
  return(err)
 }
 for (i in 1:n_trophic_group){
  checktest= (length(waste_mass[[i]])!=n_mass_list[i])
  if (checktest){
   err="inconsistent inputs of global parameters and waste mass"
   return(err)
  }
 }

 plant_biomass_series=foodweb_inputs$plant_biomass_series
 n_plant_biomass_series=length(plant_biomass_series)

 plant_senescence_parameters=foodweb_inputs$plant_senescence_parameters
 checktest= (length(plant_senescence_parameters)!=2)
 if (checktest){
  err="inconsistent input of plant senescence parameters"
  return(err)
 }
 plant_senescence_rate=plant_senescence_parameters[1]
 mass_plant_litter=plant_senescence_parameters[2]
 translation_plant_litter=Rspot(mass_plant_litter,detritus_mass,n_detritus)

 translation_juvenile_list=vector("list",n_trophic_group)
 translation_waste_list=vector("list",n_trophic_group)
 translation_carrion_list=vector("list",n_trophic_group)
 for (i in 1:n_trophic_group){
  translation_juvenile_list[[i]]=array(0,n_mass_list[i])
  translation_waste_list[[i]]=array(0,n_mass_list[i])
  translation_carrion_list[[i]]=array(0,n_mass_list[i])
  for (j in 1:n_mass_list[i]){
    translation_juvenile_list[[i]][j]=Rspot(juvenile_mass_list[[i]][j],exp(animal_logbodymass_list[[i]]),n_mass_list[i])
    translation_waste_list[[i]][j]=Rspot(waste_mass[[i]][j],detritus_mass,n_detritus)
    translation_carrion_list[[i]][j]=Rspot(animal_bodymass_average_list[[i]][j],detritus_mass,n_detritus) # carrion_mass replaced by bodymass_average
  }
 }

 pos_density_matrix=matrix(0,cols_predation_matrix,2)
 pos_predation_list=vector("list",n_trophic_group)
 pos=0
 for (i in 1:n_trophic_group){
  pos_predation_list[[i]]=array(0,n_mass_list[i])
  for (j in 1:n_mass_list[i]){
   pos_predation_list[[i]][j]=pos+j
   pos_density_matrix[(pos+j),1]=i
   pos_density_matrix[(pos+j),2]=j
  }
  pos=pos+n_mass_list[i]
 }


##########################################################
## Computing matrices of attack rates and handling times
## that are independent of population densities
## and thus only needs to be computed once
##########################################################

 attack_rate_matrix=Rcompute_attack_rate(animal_bodymass_average_list,predation_parameters,n_trophic_group,n_carnivore,n_mass_list,pos_carnivore)
 handling_time_matrix=Rcompute_handling_time(animal_bodymass_average_list,predation_parameters,n_trophic_group,n_carnivore,n_mass_list,pos_carnivore)
 herbivory_rate_list=Rcompute_herbivory_rate(animal_bodymass_average_list,pos_herbivore,herbivory_parameters,n_herbivore,n_mass_list)
 herbivory_K_list=Rcompute_herbivory_K(animal_bodymass_average_list,pos_herbivore,herbivory_parameters,n_herbivore,n_mass_list)
 detritivory_rate=Rcompute_detritivory_rates(animal_bodymass_average_list,pos_detritivore,detritivory_parameters,n_detritivore,n_mass_list)
 detritivory_attack_rate=detritivory_rate[[1]]
 detritivory_handling_time=detritivory_rate[[2]]


#########################
## FOOD WEB DYNAMICS
## predation matrix F_ij
#########################
 tot_output=n_step_dynamics/n_step_output
 i_output=1
 ad_output=list()
 det_output=list()
 output=list()
 for (istep in 1:n_step_dynamics){
  ## put waste and carrion to zero
  waste_list=Rnullify(n_trophic_group,n_mass_list)
  carrion_list=Rnullify(n_trophic_group,n_mass_list)
  food_assimilated_list=Rnullify(n_trophic_group,n_mass_list)
  food_assimilated_matrix3d=Rnullify2(n_detritivore,n_mass_list[pos_detritivore],n_detritus)

  ## read the plant (green) biomass of the corresponding time step.
  plant_biomass=plant_biomass_series[(1+((istep-1)%%n_plant_biomass_series))]

  ## herbivory and waste computation
  temp=Rherbivory_with_waste(plant_biomass,animal_density_list,pos_herbivore,herbivory_rate_list,herbivory_K_list,food_assimilated_list,animal_growth_parameters,deltaT,n_herbivore,n_mass_list,waste_list)
  food_assimilated_list=temp[[1]]
  waste_list=temp[[2]]
  plant_biomass=temp[[3]]

  if (istep==1){
   output[[1]]=food_assimilated_list
   output[[2]]=waste_list
   output[[3]]=plant_biomass
  }

  ## detritivory and waste computation
  temp=Rdetritivory_with_waste(detritus,animal_density_list,pos_detritivore,detritivory_attack_rate,detritivory_handling_time,food_assimilated_matrix3d,animal_growth_parameters,deltaT,n_detritivore,n_detritus,n_mass_list,waste_list)
  food_assimilated_matrix3d=temp[[1]]
  waste_list=temp[[2]]
  detritus=temp[[3]]
  food_assimilated_list=Rcompile_detritivore_food_eaten(food_assimilated_list,food_assimilated_matrix3d,n_detritivore,n_detritus,pos_detritivore,n_mass_list)

  if (istep==1){
   output[[4]]=food_assimilated_list
   output[[5]]=waste_list
   output[[6]]=detritus
  }

  ## computation of predation matrix F_ij
  temp=Rcomputation_predation_matrix(cols_predation_matrix,animal_density_list,attack_rate_matrix,handling_time_matrix,n_carnivore,n_mass_list,pos_density_matrix,pos_carnivore)
  predation_matrix=temp[[1]]
  sum_predation_array=temp[[2]]

  if (istep==1){
   output[[7]]=predation_matrix
   output[[8]]=sum_predation_array
  }

  ## actual predation and food eaten by each size class
  temp=Ranimal_predation_with_waste(waste_list,animal_density_list,animal_bodymass_average_list,predation_matrix,sum_predation_array,deltaT,n_trophic_group,n_carnivore,n_mass_list,animal_growth_parameters,pos_predation_list,pos_carnivore)
  animal_density_list=temp[[2]]
  food_assimilated_matrix=temp[[3]]
  waste_list=temp[[4]]
  food_assimilated_list=Rcompile_predator_food_eaten(food_assimilated_matrix,food_assimilated_list,cols_predation_matrix,pos_predation_list,n_carnivore,pos_carnivore,n_mass_list)

  if (istep==1){
   output[[9]]=animal_density_list
   output[[10]]=waste_list
   output[[11]]=food_assimilated_list
  }

  ## growth following food assimilation
  animal_density_list=Ranimal_growth_logbodymass_matrix_without_waste(animal_density_list,animal_logbodymass_list,animal_bodymass_average_list,food_assimilated_list,n_mass_list,metabolic_parameters,n_trophic_group,deltaT)

  if (istep==1){
   output[[12]]=animal_density_list
  }

  ## natural mortality
  temp=Ranimal_natural_mortality_matrix(animal_density_list,animal_bodymass_average_list,n_mass_list,mortality_parameters,deltaT,n_trophic_group)
  animal_density_list=temp[[1]]
  carrion_list=temp[[2]]

  if (istep==1){
   output[[13]]=animal_density_list
   output[[14]]=carrion_list
  }

  ## reproduction
  animal_density_list=Ranimal_reproduction_matrix(animal_density_list,animal_bodymass_average_list,juvenile_mass_list,translation_juvenile_list,reproduction_parameters,n_mass_list,deltaT,n_trophic_group)

  if (istep==1){
   output[[15]]=animal_density_list
  }

  ## update of the detritus
  detritus=Rcollect_detritus(detritus,waste_list,translation_waste_list,carrion_list,translation_carrion_list,n_trophic_group,n_mass_list)
  detritus[translation_plant_litter]=(detritus[translation_plant_litter]+plant_biomass[1]*plant_senescence_rate) ## add plant litter

  if (istep==1){
   output[[16]]=detritus
  }

  if ((istep%%n_step_output)==0){
   ## OUTPUT
   ad_output[[i_output]]=animal_density_list
   det_output[[i_output]]=detritus
   i_output=i_output+1
  }
 }

 ## OUTPUT
 list("animal_density_timeseries"=ad_output,"detritus_timeseries"=det_output,"attack_rate_matrix"=attack_rate_matrix,"handling_time_matrix"=handling_time_matrix,"herbivory_rate_list"=herbivory_rate_list,"herbivory_K_list"=herbivory_K_list,"detritivory_attack_rate"=detritivory_attack_rate,"detritivory_handling_time"=detritivory_handling_time,"output"=output,"animal_bodymass_average"=animal_bodymass_average_list,"metabolic_parameters"=metabolic_parameters)
}

