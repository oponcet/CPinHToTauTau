import os
import sys

from columnflow.util import maybe_import

np = maybe_import("numpy")
ak = maybe_import("awkward")
coffea = maybe_import("coffea")


class PhiCPNPMethod:
    def __init__(self, 
                 p4hcand1_pi : ak.Array,
                 p4hcand1_pi0: ak.Array,
                 p4hcand2_pi : ak.Array,
                 p4hcand2_pi0: ak.Array):
        """
        A class to estimate the PhiCP in decay plane method
        It will be valid only for rho and a1 decays
        arguments:
          pions and pizeros from hcand1
          pions and pizeros from hcand2
        Remarks:
          Plus and Minus signify nothing. 
          These are only the two legs of hcand.
        """
        super(PhiCPNPMethod, self).__init__()
        self.p4hcand1_pi  = p4hcand1_pi
        self.p4hcand1_pi0 = p4hcand1_pi0
        self.p4hcand2_pi  = p4hcand2_pi
        self.p4hcand2_pi0 = p4hcand2_pi0


    def get_pi_a1(self, p4_pi):
        #from IPython import embed; embed()
        p4_pi = 1 * p4_pi
        Minv1 = (p4_pi[:,:1] + p4_pi[:,1:2]).mass
        Minv2 = (p4_pi[:,:1] + p4_pi[:,2:3]).mass
        Pi = ak.where(np.abs(0.77526-Minv1) < np.abs(0.77526-Minv2), p4_pi[:,1:2], p4_pi[:,2:3])
        return Pi


    def getPhiCP_NP(self, PiPlus, PiMinus, PiZeroPlus, PiZeroMinus):
        PiPlus  = 1 * PiPlus
        PiMinus = 1 * PiMinus
        PiZeroPlus  = 1 * PiZeroPlus
        PiZeroMinus = 1 * PiZeroMinus

        #from IPython import embed; embed()
        #1/0

        ZMF = PiPlus + PiMinus
        boostv = ZMF.boostvec
        # minus side
        PiMinus_ZMF = PiMinus.boost(boostv.negative())
        PiZeroMinus_ZMF = PiZeroMinus.boost(boostv.negative())
        vecPiMinus = PiMinus_ZMF.pvec.unit
        vecPiZeroMinus = PiZeroMinus_ZMF.pvec.unit
        vecPiZeroMinustransv = (vecPiZeroMinus - vecPiMinus*(vecPiMinus.dot(vecPiZeroMinus))).unit
        # plus side
        PiPlus_ZMF = PiPlus.boost(boostv.negative())
        PiZeroPlus_ZMF = PiZeroPlus.boost(boostv.negative())
        vecPiPlus = PiPlus_ZMF.pvec.unit
        vecPiZeroPlus = PiZeroPlus_ZMF.pvec.unit
        vecPiZeroPlustransv = (vecPiZeroPlus - vecPiPlus*(vecPiPlus.dot(vecPiZeroPlus))).unit
        # Y variable
        Y1 = (PiMinus.energy - PiZeroMinus.energy)/(PiMinus.energy + PiZeroMinus.energy)
        Y2 = (PiPlus.energy - PiZeroPlus.energy)/(PiPlus.energy + PiZeroPlus.energy)
        Y  = Y1*Y2
        # angle
        acop_DP_1 = np.arccos(vecPiZeroPlustransv.dot(vecPiZeroMinustransv))
        sign_DP   = vecPiMinus.dot(vecPiZeroPlustransv.cross(vecPiZeroMinustransv))
        sign_mask = sign_DP < 0.0 
        acop_DP_2 = ak.where(sign_mask, 2*np.pi - acop_DP_1, acop_DP_1)
        Y_mask    = Y < 0.0
        acop_DP_3 = ak.where(Y_mask, acop_DP_2 + np.pi, acop_DP_2)
        mask      = Y_mask & (acop_DP_3 > 2*np.pi)
        acop_DP   = ak.where(mask, acop_DP_3 - 2*np.pi, acop_DP_3)
        
        return acop_DP
        
    
    def comp_PhiCPNP(self, tag, chbool):
        """
        Acoplanarity angle in decay plane method
        """
        PiMinus = None
        PiZeroMinus = None
        PiPlus = None
        PiZeroPlus = None
        
        p4hcand1_pi  = ak.where(chbool, self.p4hcand1_pi,  self.p4hcand1_pi[:,:0])
        p4hcand1_pi0 = ak.where(chbool, self.p4hcand1_pi0, self.p4hcand1_pi0[:,:0])
        p4hcand2_pi  = ak.where(chbool, self.p4hcand2_pi,  self.p4hcand2_pi[:,:0])
        p4hcand2_pi0 = ak.where(chbool, self.p4hcand2_pi0, self.p4hcand2_pi0[:,:0])

        if tag == "rhorho":
            PiMinus     = p4hcand1_pi
            PiZeroMinus = p4hcand1_pi0
            PiPlus      = p4hcand2_pi
            PiZeroPlus  = p4hcand2_pi0

            #from IPython import embed; embed()
            #1/0
            
        elif tag == "a1rho":
            p4hcand1_pi_3 = ak.where(ak.num(p4hcand1_pi.pt, axis=1) == 3, p4hcand1_pi, p4hcand1_pi[:,:0])
            p4hcand1_pi_1 = ak.where(ak.num(p4hcand1_pi.pt, axis=1) == 1, p4hcand1_pi, p4hcand1_pi[:,:0])
            PiMinus = ak.where(ak.num(p4hcand1_pi.pt, axis=1) == 3,
                               self.get_pi_a1(p4hcand1_pi_3),
                               p4hcand1_pi_1)
            PiZeroMinus = ak.where(ak.num(p4hcand1_pi.pt, axis=1) == 3, 
                                   p4hcand1_pi[:,:1],
                                   p4hcand1_pi0)
            p4hcand2_pi_3 = ak.where(ak.num(p4hcand2_pi.pt, axis=1) == 3, p4hcand2_pi, p4hcand2_pi[:,:0])
            p4hcand2_pi_1 = ak.where(ak.num(p4hcand2_pi.pt, axis=1) == 1, p4hcand2_pi, p4hcand2_pi[:,:0])
            PiPlus = ak.where(ak.num(p4hcand2_pi.pt, axis=1) == 3,
                              self.get_pi_a1(p4hcand2_pi_3),
                              p4hcand2_pi_1)
            PiZeroPlus = ak.where(ak.num(p4hcand2_pi.pt, axis=1) == 3, 
                                  p4hcand2_pi[:,:1],
                                  p4hcand2_pi0)

            #from IPython import embed; embed()
            #1/0
            
        elif tag == "a1a1":
            PiMinus     = self.get_pi_a1(p4hcand1_pi)
            PiZeroMinus = p4hcand1_pi[:,:1]
            PiPlus      = self.get_pi_a1(p4hcand2_pi)
            PiZeroPlus  = p4hcand2_pi[:,:1]

            #from IPython import embed; embed()
            #1/0
            
        else:
            raise RuntimeWarning (f"{self.tag} is not suitable for Neutral Pion method")

        phicp = self.getPhiCP_NP(PiPlus, PiMinus, PiZeroPlus, PiZeroMinus)
        
        return phicp
