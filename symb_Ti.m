Ti = [
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               -Lstance*sin(gamma1)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                Lstance*cos(gamma1)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               -Lstance*sin(gamma1)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                Lstance*cos(gamma1)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             Lhip/2
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               -Lstance*sin(gamma1)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                Lstance*cos(gamma1)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               Lhip
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    - Lstance*sin(gamma1) - cgThigh*cos(alpha2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      Lstance*cos(gamma1) - cgThigh*cos(alpha2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         Lhip + cgThigh*sin(alpha2)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     - Lstance*sin(gamma1) - Lthigh*cos(alpha2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Lstance*cos(gamma1) - Lthigh*cos(alpha2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          Lhip + Lthigh*sin(alpha2)
                                                                                                                                                                                                                                                   cgShank*(sin(gamma3)*(cos(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2)) + sin(alpha2)*sin(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - cos(alpha2)*cos(gamma3)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - Lstance*sin(gamma1) - Lthigh*cos(alpha2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))
                                                                                                                                                                                                                                                   Lstance*cos(gamma1) - cgShank*(sin(gamma3)*(cos(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1)) - sin(alpha2)*sin(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) + cos(alpha2)*cos(gamma3)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) - Lthigh*cos(alpha2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Lhip + Lthigh*sin(alpha2) + cgShank*(cos(gamma3)*sin(alpha2) + cos(alpha2)*sin(beta2)*sin(gamma3))
                                                                                                                                                                                                                                                    Lshank*(sin(gamma3)*(cos(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2)) + sin(alpha2)*sin(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - cos(alpha2)*cos(gamma3)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - Lstance*sin(gamma1) - Lthigh*cos(alpha2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))
                                                                                                                                                                                                                                                    Lstance*cos(gamma1) - Lshank*(sin(gamma3)*(cos(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1)) - sin(alpha2)*sin(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) + cos(alpha2)*cos(gamma3)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) - Lthigh*cos(alpha2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Lhip + Lthigh*sin(alpha2) + Lshank*(cos(gamma3)*sin(alpha2) + cos(alpha2)*sin(beta2)*sin(gamma3))
 cgFoot*(sin(gamma4)*(cos(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2)) + sin(alpha2)*sin(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - cos(alpha2)*cos(gamma4)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - Lstance*sin(gamma1) + Lshank*(sin(gamma3)*(cos(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2)) + sin(alpha2)*sin(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - cos(alpha2)*cos(gamma3)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - Lthigh*cos(alpha2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))
 Lstance*cos(gamma1) - cgFoot*(sin(gamma4)*(cos(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1)) - sin(alpha2)*sin(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) + cos(alpha2)*cos(gamma4)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) - Lshank*(sin(gamma3)*(cos(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1)) - sin(alpha2)*sin(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) + cos(alpha2)*cos(gamma3)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) - Lthigh*cos(alpha2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))
                                                                                                                                                                                                                                                                                                                                                                                                                          Lhip + Lthigh*sin(alpha2) + Lshank*(cos(gamma3)*sin(alpha2) + cos(alpha2)*sin(beta2)*sin(gamma3)) + cgFoot*(cos(gamma4)*sin(alpha2) + cos(alpha2)*sin(beta2)*sin(gamma4))
  Lfoot*(sin(gamma4)*(cos(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2)) + sin(alpha2)*sin(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - cos(alpha2)*cos(gamma4)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - Lstance*sin(gamma1) + Lshank*(sin(gamma3)*(cos(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2)) + sin(alpha2)*sin(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - cos(alpha2)*cos(gamma3)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))) - Lthigh*cos(alpha2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1))
  Lstance*cos(gamma1) - Lfoot*(sin(gamma4)*(cos(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1)) - sin(alpha2)*sin(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) + cos(alpha2)*cos(gamma4)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) - Lshank*(sin(gamma3)*(cos(beta2)*(cos(gamma1)*sin(gamma2) + cos(gamma2)*sin(gamma1)) - sin(alpha2)*sin(beta2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) + cos(alpha2)*cos(gamma3)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))) - Lthigh*cos(alpha2)*(sin(gamma1)*sin(gamma2) - cos(gamma1)*cos(gamma2))
                                                                                                                                                                                                                                                                                                                                                                                                                           Lhip + Lthigh*sin(alpha2) + Lfoot*(cos(gamma4)*sin(alpha2) + cos(alpha2)*sin(beta2)*sin(gamma4)) + Lshank*(cos(gamma3)*sin(alpha2) + cos(alpha2)*sin(beta2)*sin(gamma3))
 
];