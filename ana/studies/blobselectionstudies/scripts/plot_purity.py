import ROOT

ROOT.gROOT.SetBatch()

f=ROOT.TFile("/minerva/data/users/tejinc/CCQENu/antinu/hists/studies/recoil_study_minervame5A.root")

rb=10

h_num = f.Get("h_blobTotalE_parentpid2112_pid2212_qelike")
h_num.Rebin(rb)
h_num.GetXaxis().SetRangeUser(0,300)
h_num.GetXaxis().SetTitle("Candidate Energy (MeV)")
h_denom = f.Get("h_blobTotalE_total_qelike")
h_denom.Rebin(rb)

print( h_num.Integral()/h_denom.Integral()*100 )

h_num.Divide( h_num, h_denom )

c=ROOT.TCanvas()

h_num.Draw("hist")

c.Print("eff_neutron_sel.pdf")

