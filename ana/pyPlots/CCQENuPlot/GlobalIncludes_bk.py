import ROOT
def enum(*args):
    enums = dict(zip(args, range(len(args))))
    return type('Enum', (), enums)

kHistos = enum(
    "kData",
    "kMC",
    "kQE",
    "kQE_H",
    "kQE_OTH",
    "kQENot",
    "k2p2h",
    "kRES",
    "kDIS",
    "kOTH",
    "kQENot_PositivePions",
    "kQENot_PositivePions_TrueMichel",
    "kQENot_PositivePions_TrueMichelNot",
    "kQENot_NegativePions",
    "kQENot_NeutralPions",
    "kQENot_NoPions",
    "kQELike_QE",
    "kQELike_QE_H",
    "kQELike_QE_OTH",
    "kQELike_2p2h",
    "kQELike_RES",
    "kQELike_DIS",
    "kQELike_OTH",
    "kQELike",
    "kQELikeNot",
    "kQELikeNot_NoPions",
    "kQELikeNot_PositivePions",
    "kQELikeNot_PositivePions_TrueMichel",
    "kQELikeNot_PositivePions_TrueMichelNot",
    "kQELikeNot_NegativePions",
    "kQELikeNot_NeutralPions",
    "kQELikeNot_ChargedPions",
    "kQELikeNot_ChargedPions_TrueMichel",
    "kQELikeNot_ChargedPions_TrueMichelNot",
    "kQELikeNot_ChargedPionsNot",
    "kQELike_QENot",
    "kQELikeNot_QE",
    "kQELikeNot_QENot",
    "kQELikeNot_2p2h",
    "kQELikeNot_SinglePion",
    "kQELikeNot_SingleChargedPion",
    "kQELikeNot_SingleNeutralPion",
    "kQELikeNot_MultiPion",
    "kQELikeNot_Other",
    "kQELikeNot_DIS",
    "kQELikeNot_RES",
    "kQELikeNot_Coh",
    "kQELike_RES_NeutronFSI",
    "kQELike_RES_ProtonFSI",
    "kQELike_RES_NeutronIS",
    "kQELike_RES_ProtonIS",
    "kQELike_QE_ProtonFSI",
    "kQELike_QE_NeutronFSI",
    "kQELike_2p2h_np",
    "kQELike_2p2h_nn",
    "kBlob_Nuc_QELike_QE_H",
    "kBlob_Nuc_QELike_QE_OTH",
    "kBlob_Nuc_QELike_RES",
    "kBlob_Nuc_QELike_2p2h",
    "kBlob_Nuc_QELike_DIS",
    "kBlob_Nuc_QELikeNot",
    "kBlob_EM_QELike_QE_H",
    "kBlob_EM_QELike_QE_OTH",
    "kBlob_EM_QELike_RES",
    "kBlob_EM_QELike_2p2h",
    "kBlob_EM_QELike_DIS",
    "kBlob_EM_QELikeNot",
    "kBlob_OTH_QELike_QE_H",
    "kBlob_OTH_QELike_QE_OTH",
    "kBlob_OTH_QELike_RES",
    "kBlob_OTH_QELike_2p2h",
    "kBlob_OTH_QELike_DIS",
    "kBlob_OTH_QELikeNot",
    "kBlob_None")

names=(
    "data",
    "mc",
    "qe",
    "qe_h",
    "qe_oth",
    "qenot",
    "2p2h",
    "res",
    "dis",
    "oth",
    "qenot_positivepions",
    "qenot_positivepions_truemichel",
    "qenot_positivepions_truemichelnot",
    "qenot_negativepions",
    "qenot_neutralpions",
    "qenot_nopions",
    "qelike_qe",
    "qelike_qe_h",
    "qelike_qe_oth",
    "qelike_2p2h",
    "qelike_res",
    "qelike_dis",
    "qelike_oth",
    "qelike",
    "qelikenot",
    "qelikenot_nopions",
    "qelikenot_positivepions",
    "qelikenot_positivepions_truemichel",
    "qelikenot_positivepions_truemichelnot",
    "qelikenot_negativepions",
    "qelikenot_neutralpions",
    "qelikenot_chargedpions",
    "qelikenot_chargedpions_truemichel",
    "qelikenot_chargedpions_truemichelnot",
    "qelikenot_chargedpionsnot",
    "qelike_qenot",
    "qelikenot_qe",
    "qelikenot_qenot",
    "qelikenot_2p2h",
    "qelikenot_singlepion",
    "qelikenot_singlechargedpion",
    "qelikenot_singleneutralpion",
    "qelikenot_multipion",
    "qelikenot_not_scp_snp_mp",
    "qelikenot_dis",
    "qelikenot_res",
    "qelikenot_coh",
    "qelike_res_neutron_fsi",
    "qelike_res_proton_fsi",
    "qelike_res_neutron_is",
    "qelike_res_proton_is",
    "qelike_qe_proton_fsi",
    "qelike_qe_neutron_fsi",
    "qelike_2p2h_np",
    "qelike_2p2h_nn",


    "blob_nuc_qelike_qe_h",
    "blob_nuc_qelike_qe_oth",
    "blob_nuc_qelike_res",
    "blob_nuc_qelike_2p2h",
    "blob_nuc_qelike_dis",
    "blob_nuc_qelikenot",
    "blob_em_qelike_qe_h",
    "blob_em_qelike_qe_oth",
    "blob_em_qelike_res",
    "blob_em_qelike_2p2h",
    "blob_em_qelike_dis",
    "blob_em_qelikenot",
    "blob_oth_qelike_qe_h",
    "blob_oth_qelike_qe_oth",
    "blob_oth_qelike_res",
    "blob_oth_qelike_2p2h",
    "blob_oth_qelike_dis",
    "blob_oth_qelikenot",
    "blob_none" )

colors = (
      ROOT.kBlack,      #Data
      ROOT.kRed,        #MC
      ROOT.kBlue-7,     #QE
      ROOT.kBlue-5,     #QE_H
      ROOT.kBlue-2,     #QE_OTH
      ROOT.kGreen-6,    #QENot
      ROOT.kViolet-8,   #2p2h
      ROOT.kOrange+1,   #RES
      ROOT.kGreen-2,    #DIS
      ROOT.kRed+4,      #OTH
      ROOT.kBlue-7,     #QENot_PositivePions
      ROOT.kOrange-3,   #QENot_PositivePions_TrueMichel
      ROOT.kMagenta+3,  #QENot_PositivePions_TrueMichelNot
      ROOT.kGreen-3,    #QENot_NegativePions
      ROOT.kRed-4,      #QENot_pi0
      ROOT.kRed+4,      #QENot_nopions
      ROOT.kGreen+1,    #QELike_QE
      ROOT.kRed,        #QELike_QE_H
      ROOT.kGreen+1,    #QELike_QE_OTH
      ROOT.kPink-9,     #QELike_2p2h
      ROOT.kAzure+1,    #QELike_RES
      ROOT.kOrange,     #QELike_DIS
      ROOT.kBlack,      #QELike_OTH
      ROOT.kMagenta-4,  #QELike
      ROOT.kBlue-7,     #QELikeNot
      ROOT.kRed+4,      #QELikeNot_NoPions
      ROOT.kBlue-7,     #QELikeNot_PositivePions
      ROOT.kOrange-3,   #QELikeNot_PositivePions_TrueMichel
      ROOT.kMagenta+3,  #QELikeNot_PositivePions_TrueMichelNot
      ROOT.kGreen-3,    #QELikeNot_NegativePions
      ROOT.kRed-4,      #QELikeNot_pi0
      ROOT.kBlue-7,     #QELikeNot_ChargedPions
      ROOT.kOrange-3,   #QELikeNot_ChargedPions_TrueMichel
      ROOT.kMagenta+3,  #QELikeNot_ChargedPions_TrueMichelNot
      ROOT.kRed-4,      #QELikeNot_ChargedPionsNot
      ROOT.kYellow-3,   #QELike_QENot
      ROOT.kTeal+5,     #QELikeNot_QE
      ROOT.kRed-7,      #QELikeNot_QENot
      ROOT.kViolet+10,  #QELikeNot_2p2h
      ROOT.kCyan+10,    #QELikeNot_SinglePions
      ROOT.kBlue-7,     #QELikeNot_SingleChargedPions
      ROOT.kGreen-3,    #QELikeNot_SingleNeutralPions
      ROOT.kRed-4,      #QELikeNot_MultiPion
      ROOT.kGreen+2,    #QELikeNot_Other
      ROOT.kOrange-9,   #QELikeNot_DIS
      ROOT.kGreen+4,    #QELikeNot_RES
      ROOT.kGreen-2    #QELikeNot_COH
    )















