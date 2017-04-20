/* These Efficiencies have been derived for Summer15ttbar events and should 
   be used only for the same MC samples or for events with similar topology */  


float BTagSFUtil::TagEfficiencyB(float JetPt, float JetEta) {
  
  if (TaggerOP=="CSVv2M") {
    if (JetPt > 20 && JetPt <= 30){
      if      (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.542647;
      else if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.526381;
      else if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.46804;
      else if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.370068;
    }

    else if (JetPt > 30 && JetPt <= 50){
      if      (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.643642;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.630665;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.584715;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.498335;
    }

    else if (JetPt > 50 && JetPt <= 70){
      if      (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.678361;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.666078;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.615402;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.528537;
    }

    else if (JetPt > 70 && JetPt <= 100){
      if      (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.691631;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.683235;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.627146;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.533232;
    }

    else if (JetPt > 100 && JetPt <= 140){
      if      (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.691601;
      else if (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.686815;
      else if (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.61803;
      else if (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.523252;
    }

    else if (JetPt > 140 && JetPt <= 200){
      if      (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.675449;
      else if (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.673037;
      else if (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.602255;
      else if (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.513437;
    }

    else if (JetPt > 200 && JetPt <= 300){
      if      (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.63579;
      else if (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.632677;
      else if (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.54259;
      else if (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.456817;
    }

    else if (JetPt > 300 && JetPt <= 3000){
      if      (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.547882;
      else if (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.561086;
      else if (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.496071;
      else if (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.446656;
    }

  }


  
  if (TaggerOP=="CSVv2T") { 
    if (JetPt > 20.00000 && JetPt <= 40.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.59234 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.57961;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.52772;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.43535;  
    } 
    else if (JetPt > 40.00000 && JetPt <= 60.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.66022 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.64971;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.60206;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.52218;  
    } 
    else if (JetPt > 60.00000 && JetPt <= 80.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68184 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.67571;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.62020;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.53328;  
    } 
    else if (JetPt > 80.00000 && JetPt <= 100.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68998 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.68509;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.63066;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.53861;  
    } 
    else if (JetPt > 100.00000 && JetPt <= 120.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.68730 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.68463;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.61622;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.52498;  
    } 
    else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.66520 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.66545;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.59499;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.51355;  
    } 
  } 
  return 1.;
}

float BTagSFUtil::TagEfficiencyC(float JetPt, float JetEta) {
  
  if (TaggerOP=="CSVv2M") {
    if (JetPt > 20 && JetPt <= 30){
      if      (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.116727;
      else if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.109121;
      else if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0984467;
      else if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0732115;
    }

    else if (JetPt > 30 && JetPt <= 50){
      if      (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.130428;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.124702;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.12614;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.103106;
    }

    else if (JetPt > 50 && JetPt <= 70){
      if      (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.126891;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.121356;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.122247;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0991674;
    }

    else if (JetPt > 70 && JetPt <= 100){
      if      (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.131153;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.129504;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.127816;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.104056;
    }

    else if (JetPt > 100 && JetPt <= 140){
      if      (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.137708;
      else if (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.140439;
      else if (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.129276;
      else if (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.110438;
    }

    else if (JetPt > 140 && JetPt <= 200){
      if      (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.135622;
      else if (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.143716;
      else if (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.135302;
      else if (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.122262;
    }

    else if (JetPt > 200 && JetPt <= 300){
      if      (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.127932;
      else if (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.140633;
      else if (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.120885;
      else if (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.113195;
    }

    else if (JetPt > 300 && JetPt <= 3000){
      if      (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.11432;
      else if (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.138309;
      else if (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.135161;
      else if (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.13728;
    }

  }

  
  if (TaggerOP=="CSVv2T") { 
    if (JetPt > 20.00000 && JetPt <= 40.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12617 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11974;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11290;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.08717;  
    } 
    else if (JetPt > 40.00000 && JetPt <= 60.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12341 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.11844;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.11937;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.09907;  
    } 
    else if (JetPt > 60.00000 && JetPt <= 80.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.12777 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.12588;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12337;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10186;  
    } 
    else if (JetPt > 80.00000 && JetPt <= 100.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.13114 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.13176;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12900;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10778;  
    } 
    else if (JetPt > 100.00000 && JetPt <= 120.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.13586 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.13863;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.12548;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.10638;  
    } 
    else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.13414 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.14295;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.13146;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.12254;  
    } 
  } 
  return 1.;
}
float BTagSFUtil::TagEfficiencyLight(float JetPt, float JetEta) {
  

  if (TaggerOP=="CSVv2M") {
    if (JetPt > 20 && JetPt <= 30){
      if      (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.00989268;
      else if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0110037;
      else if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0120286;
      else if (fabs(JetEta) > 0 && fabs(JetEta) <= 0.6) return 0.0115793;
    }

    else if (JetPt > 30 && JetPt <= 50){
      if      (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.00945715;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0111287;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0147667;
      else if (fabs(JetEta) > 0.6 && fabs(JetEta) <= 1.2) return 0.0174294;
    }

    else if (JetPt > 50 && JetPt <= 70){
      if      (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.00784805;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.00935706;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0117512;
      else if (fabs(JetEta) > 1.2 && fabs(JetEta) <= 1.8) return 0.0148856;
    }

    else if (JetPt > 70 && JetPt <= 100){
      if      (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.00806538;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.00978417;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0119189;
      else if (fabs(JetEta) > 1.8 && fabs(JetEta) <= 2.4) return 0.0150979;
    }

    else if (JetPt > 100 && JetPt <= 140){
      if      (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.00820627;
      else if (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.0101786;
      else if (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.0124413;
      else if (fabs(JetEta) > 2.4 && fabs(JetEta) <= 3) return 0.0171922;
    }

    else if (JetPt > 140 && JetPt <= 200){
      if      (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.00869274;
      else if (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.011655;
      else if (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.0166625;
      else if (fabs(JetEta) > 3 && fabs(JetEta) <= 3.6) return 0.0243483;
    }

    else if (JetPt > 200 && JetPt <= 300){
      if      (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.00982392;
      else if (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.0140469;
      else if (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.0191645;
      else if (fabs(JetEta) > 3.6 && fabs(JetEta) <= 4.2) return 0.0257676;
    }

    else if (JetPt > 300 && JetPt <= 3000){
      if      (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.0140646;
      else if (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.0215447;
      else if (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.0326195;
      else if (fabs(JetEta) > 4.2 && fabs(JetEta) <= 4.8) return 0.0362126;
    }

  }

  
  if (TaggerOP=="CSVv2T") { 
    if (JetPt > 20.00000 && JetPt <= 40.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00953 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01103;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01286;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01357;  
    } 
    else if (JetPt > 40.00000 && JetPt <= 60.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00690 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00841;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01088;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01451;  
    } 
    else if (JetPt > 60.00000 && JetPt <= 80.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00656 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00807;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00992;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01298;  
  } 
    else if (JetPt > 80.00000 && JetPt <= 100.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00685 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00847;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01069;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01411;  
    } 
    else if (JetPt > 100.00000 && JetPt <= 120.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00675 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.00837;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.00987;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.01367;  
    } 
    else if (JetPt > 120.00000 && JetPt <= 3000.00000){  
      if      (fabs(JetEta) > 0.00000   && fabs(JetEta) <= 0.60000) return 0.00816 ;  
      else if (fabs(JetEta) > 0.60000 && fabs(JetEta) <= 1.20000) return 0.01119;  
      else if (fabs(JetEta) > 1.20000 && fabs(JetEta) <= 1.80000) return 0.01543;  
      else if (fabs(JetEta) > 1.80000 && fabs(JetEta) <= 2.40000) return 0.02225;  
    } 
  } 
  return 1.;
}
