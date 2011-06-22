// Made on Mon, Jun 20, 2011 at 9:00 PM using
//   python estimatePileupD.py --maxPileupBin=50 -i ../DataMCSpectraComparison/ana_datamc_Run2011AMuonsOnly/ana_datamc_data.forlumi.json pudist.root --debugLumi
// for 746/pb: see python/goodlumis.py in CVS from that date.

{
//=========Macro generated from canvas: Canvas 1/Canvas 1
//=========  (Wed Jun 22 07:11:45 2011) by ROOT version5.27/04
   TCanvas *Canvas 1 = new TCanvas("Canvas 1", "Canvas 1",735,108,1137,939);
   Canvas->Range(-6.875,-1.495759e+007,56.875,1.346183e+008);
   Canvas->SetBorderSize(2);
   Canvas->SetFrameFillColor(0);
   
   TH1D *pileup = new TH1D("pileup","pileup",51,-0.5,50.5);
   pileup->SetBinContent(1,1.070966e+007);
   pileup->SetBinContent(2,2.078337e+007);
   pileup->SetBinContent(3,4.838647e+007);
   pileup->SetBinContent(4,8.042218e+007);
   pileup->SetBinContent(5,1.049422e+008);
   pileup->SetBinContent(6,1.139626e+008);
   pileup->SetBinContent(7,1.069211e+008);
   pileup->SetBinContent(8,8.895942e+007);
   pileup->SetBinContent(9,6.690898e+007);
   pileup->SetBinContent(10,4.616308e+007);
   pileup->SetBinContent(11,2.955305e+007);
   pileup->SetBinContent(12,1.771738e+007);
   pileup->SetBinContent(13,1.002196e+007);
   pileup->SetBinContent(14,5382477);
   pileup->SetBinContent(15,2759227);
   pileup->SetBinContent(16,1356246);
   pileup->SetBinContent(17,641709.9);
   pileup->SetBinContent(18,293273.8);
   pileup->SetBinContent(19,129849.9);
   pileup->SetBinContent(20,55844.82);
   pileup->SetBinContent(21,23382.74);
   pileup->SetBinContent(22,9551.016);
   pileup->SetBinContent(23,3812.41);
   pileup->SetBinContent(24,1489.332);
   pileup->SetBinContent(25,570.1266);
   pileup->SetBinContent(26,214.0871);
   pileup->SetBinContent(27,78.92504);
   pileup->SetBinContent(28,28.58447);
   pileup->SetBinContent(29,10.17541);
   pileup->SetBinContent(30,3.561506);
   pileup->SetBinContent(31,1.225961);
   pileup->SetBinContent(32,0.4150875);
   pileup->SetBinContent(33,0.1382444);
   pileup->SetBinContent(34,0.04528982);
   pileup->SetBinContent(35,0.01459426);
   pileup->SetBinContent(36,0.004625555);
   pileup->SetBinContent(37,0.001441825);
   pileup->SetBinContent(38,0.0004419702);
   pileup->SetBinContent(39,0.0001332207);
   pileup->SetBinContent(40,3.948381e-005);
   pileup->SetBinContent(41,1.150564e-005);
   pileup->SetBinContent(42,3.296313e-006);
   pileup->SetBinContent(43,9.28456e-007);
   pileup->SetBinContent(44,2.571019e-007);
   pileup->SetBinContent(45,6.999444e-008);
   pileup->SetBinContent(46,1.873469e-008);
   pileup->SetBinContent(47,4.930277e-009);
   pileup->SetBinContent(48,1.275739e-009);
   pileup->SetBinContent(49,3.245995e-010);
   pileup->SetBinContent(50,2.097218e-007);
   pileup->SetEntries(4.850141e+009);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.835,0.98,0.995,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(2);
   ptstats->SetFillColor(19);
   ptstats->SetTextAlign(12);
   TText *text = ptstats->AddText("pileup");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries =  4.850141e+009");
   text = ptstats->AddText("Mean  =  5.717");
   text = ptstats->AddText("RMS   =  2.724");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   pileup->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(pileup->GetListOfFunctions());
   pileup->Draw("");
   
   TPaveText *pt = new TPaveText(0.01,0.9420032,0.112397,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(2);
   pt->SetFillColor(19);
   text = pt->AddText("pileup");
   pt->Draw();
   Canvas->Modified();
   Canvas 1->cd();
   Canvas 1->SetSelected(Canvas 1);
}
