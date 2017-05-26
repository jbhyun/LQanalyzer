float BTagSFUtil::TagEfficiencyB(float JetPt, float JetEta) {


  if (TaggerOP=="CSVv2M") {
    if (JetPt > 20 && JetPt <= 30){
      if      (fabs(JetEta) <= 0.6) return 0.542994;
      else if (fabs(JetEta) <= 1.2) return 0.526664;
      else if (fabs(JetEta) <= 1.8) return 0.468419;
      else if (fabs(JetEta) <= 2.4) return 0.370355;
    }

    else if (JetPt > 30 && JetPt <= 50){
      if      (fabs(JetEta) <= 0.6) return 0.644011;
      else if (fabs(JetEta) <= 1.2) return 0.631047;
      else if (fabs(JetEta) <= 1.8) return 0.585269;
      else if (fabs(JetEta) <= 2.4) return 0.499157;
    }

    else if (JetPt > 50 && JetPt <= 70){
      if      (fabs(JetEta) <= 0.6) return 0.678655;
      else if (fabs(JetEta) <= 1.2) return 0.666408;
      else if (fabs(JetEta) <= 1.8) return 0.615887;
      else if (fabs(JetEta) <= 2.4) return 0.529223;
    }

    else if (JetPt > 70 && JetPt <= 100){
      if      (fabs(JetEta) <= 0.6) return 0.691881;
      else if (fabs(JetEta) <= 1.2) return 0.683536;
      else if (fabs(JetEta) <= 1.8) return 0.627555;
      else if (fabs(JetEta) <= 2.4) return 0.533835;
    }

    else if (JetPt > 100 && JetPt <= 140){
      if      (fabs(JetEta) <= 0.6) return 0.691809;
      else if (fabs(JetEta) <= 1.2) return 0.687077;
      else if (fabs(JetEta) <= 1.8) return 0.618308;
      else if (fabs(JetEta) <= 2.4) return 0.523716;
    }

    else if (JetPt > 140 && JetPt <= 200){
      if      (fabs(JetEta) <= 0.6) return 0.675649;
      else if (fabs(JetEta) <= 1.2) return 0.67315;
      else if (fabs(JetEta) <= 1.8) return 0.602469;
      else if (fabs(JetEta) <= 2.4) return 0.513714;
    }

    else if (JetPt > 200 && JetPt <= 300){
      if      (fabs(JetEta) <= 0.6) return 0.635681;
      else if (fabs(JetEta) <= 1.2) return 0.632741;
      else if (fabs(JetEta) <= 1.8) return 0.542772;
      else if (fabs(JetEta) <= 2.4) return 0.45699;
    }

    else if (JetPt > 300){
      if      (fabs(JetEta) <= 0.6) return 0.547572;
      else if (fabs(JetEta) <= 1.2) return 0.56159;
      else if (fabs(JetEta) <= 1.8) return 0.495854;
      else if (fabs(JetEta) <= 2.4) return 0.445693;
    }

  }
}


float BTagSFUtil::TagEfficiencyC(float JetPt, float JetEta) {


  if (TaggerOP=="CSVv2M") {
    if (JetPt > 20 && JetPt <= 30){
      if      (fabs(JetEta) <= 0.6) return 0.116709;
      else if (fabs(JetEta) <= 1.2) return 0.109228;
      else if (fabs(JetEta) <= 1.8) return 0.098289;
      else if (fabs(JetEta) <= 2.4) return 0.0731406;
    }

    else if (JetPt > 30 && JetPt <= 50){
      if      (fabs(JetEta) <= 0.6) return 0.130515;
      else if (fabs(JetEta) <= 1.2) return 0.124753;
      else if (fabs(JetEta) <= 1.8) return 0.12634;
      else if (fabs(JetEta) <= 2.4) return 0.103325;
    }

    else if (JetPt > 50 && JetPt <= 70){
      if      (fabs(JetEta) <= 0.6) return 0.126951;
      else if (fabs(JetEta) <= 1.2) return 0.121492;
      else if (fabs(JetEta) <= 1.8) return 0.12241;
      else if (fabs(JetEta) <= 2.4) return 0.0992856;
    }

    else if (JetPt > 70 && JetPt <= 100){
      if      (fabs(JetEta) <= 0.6) return 0.131285;
      else if (fabs(JetEta) <= 1.2) return 0.129668;
      else if (fabs(JetEta) <= 1.8) return 0.127941;
      else if (fabs(JetEta) <= 2.4) return 0.10436;
    }

    else if (JetPt > 100 && JetPt <= 140){
      if      (fabs(JetEta) <= 0.6) return 0.137849;
      else if (fabs(JetEta) <= 1.2) return 0.140735;
      else if (fabs(JetEta) <= 1.8) return 0.129565;
      else if (fabs(JetEta) <= 2.4) return 0.110414;
    }

    else if (JetPt > 140 && JetPt <= 200){
      if      (fabs(JetEta) <= 0.6) return 0.135571;
      else if (fabs(JetEta) <= 1.2) return 0.144002;
      else if (fabs(JetEta) <= 1.8) return 0.135353;
      else if (fabs(JetEta) <= 2.4) return 0.122589;
    }

    else if (JetPt > 200 && JetPt <= 300){
      if      (fabs(JetEta) <= 0.6) return 0.128062;
      else if (fabs(JetEta) <= 1.2) return 0.141174;
      else if (fabs(JetEta) <= 1.8) return 0.120956;
      else if (fabs(JetEta) <= 2.4) return 0.113396;
    }

    else if (JetPt > 300){
      if      (fabs(JetEta) <= 0.6) return 0.114269;
      else if (fabs(JetEta) <= 1.2) return 0.137899;
      else if (fabs(JetEta) <= 1.8) return 0.135165;
      else if (fabs(JetEta) <= 2.4) return 0.137883;
    }

  }
}


float BTagSFUtil::TagEfficiencyLight(float JetPt, float JetEta) {


  if (TaggerOP=="CSVv2M") {
    if (JetPt > 20 && JetPt <= 30){
      if      (fabs(JetEta) <= 0.6) return 0.00990635;
      else if (fabs(JetEta) <= 1.2) return 0.0109924;
      else if (fabs(JetEta) <= 1.8) return 0.0120132;
      else if (fabs(JetEta) <= 2.4) return 0.0115162;
    }

    else if (JetPt > 30 && JetPt <= 50){
      if      (fabs(JetEta) <= 0.6) return 0.00945299;
      else if (fabs(JetEta) <= 1.2) return 0.0111265;
      else if (fabs(JetEta) <= 1.8) return 0.014771;
      else if (fabs(JetEta) <= 2.4) return 0.0174325;
    }

    else if (JetPt > 50 && JetPt <= 70){
      if      (fabs(JetEta) <= 0.6) return 0.00782955;
      else if (fabs(JetEta) <= 1.2) return 0.00934679;
      else if (fabs(JetEta) <= 1.8) return 0.0117663;
      else if (fabs(JetEta) <= 2.4) return 0.0148924;
    }

    else if (JetPt > 70 && JetPt <= 100){
      if      (fabs(JetEta) <= 0.6) return 0.00804995;
      else if (fabs(JetEta) <= 1.2) return 0.00977188;
      else if (fabs(JetEta) <= 1.8) return 0.0118909;
      else if (fabs(JetEta) <= 2.4) return 0.0150796;
    }

    else if (JetPt > 100 && JetPt <= 140){
      if      (fabs(JetEta) <= 0.6) return 0.0081993;
      else if (fabs(JetEta) <= 1.2) return 0.0101671;
      else if (fabs(JetEta) <= 1.8) return 0.0124102;
      else if (fabs(JetEta) <= 2.4) return 0.0171818;
    }

    else if (JetPt > 140 && JetPt <= 200){
      if      (fabs(JetEta) <= 0.6) return 0.00868441;
      else if (fabs(JetEta) <= 1.2) return 0.011653;
      else if (fabs(JetEta) <= 1.8) return 0.0166419;
      else if (fabs(JetEta) <= 2.4) return 0.0243221;
    }

    else if (JetPt > 200 && JetPt <= 300){
      if      (fabs(JetEta) <= 0.6) return 0.00985362;
      else if (fabs(JetEta) <= 1.2) return 0.0140147;
      else if (fabs(JetEta) <= 1.8) return 0.0191696;
      else if (fabs(JetEta) <= 2.4) return 0.0258203;
    }

    else if (JetPt > 300){
      if      (fabs(JetEta) <= 0.6) return 0.0140786;
      else if (fabs(JetEta) <= 1.2) return 0.0215001;
      else if (fabs(JetEta) <= 1.8) return 0.0327833;
      else if (fabs(JetEta) <= 2.4) return 0.0364012;
    }

  }
}
