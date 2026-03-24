import ROOT as R

def main():
    rootFile = "/gpfs/group/belle/users/wangz/data_gMC_belle1/KsSpinAlignment_v2.3.2_qqbar_svd2/svd2_reco_processed.root"
    #df = R.RDataFrame("event", "temp_bin_10.root")
    df = R.RDataFrame("event", rootFile)
    df =  df.Filter("Ks_M.size() != Ks_weight.size()", "not match size")
    df.Report().Print()

main()