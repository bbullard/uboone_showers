#include <TH1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <time.h>
#include <ctime>
#include <cstdlib>

using namespace std;
using namespace TMath;

/*
 * Returns T/F if input point is contained within TPC volume
 */ 
bool PointInTPCVolume(array<double, 3> p)
{
	double x_fid_min = 0.0;
        double x_fid_max = 256.35;
        double y_fid_min = -116.5;
        double y_fid_max = 116.5;
        double z_fid_min = 0.0;
        double z_fid_max = 1036.8;

	bool inX = (x_fid_min <= p[0]) and (p[0] <= x_fid_max);
	bool inY = (y_fid_min <= p[1]) and (p[1] <= y_fid_max);
	bool inZ = (z_fid_min <= p[2]) and (p[2] <= z_fid_max);	
	
	return (inX and inY and inZ);
}


/*
 * Takes as input an a point, a unit array, and an angle.
 * Returns the rotated point around the origin and the array
 */
array<double, 3> RotatePoint(array<double, 3> p, array<double, 3> v, double theta)
{
	array<double, 3> result = {{0.,0.,0.}};
	
	result[0] = v[0]*(v[0]*p[0] + v[1]*p[1] + v[2]*p[2])*(1-Cos(theta)) + p[0]*Cos(theta) + (-v[2]*p[1] + v[1]*p[2])*Sin(theta);
	result[1] = v[1]*(v[0]*p[0] + v[1]*p[1] + v[2]*p[2])*(1-Cos(theta)) + p[1]*Cos(theta) + (v[2]*p[0] - v[0]*p[2])*Sin(theta); 	
	result[2] = v[2]*(v[0]*p[0] + v[1]*p[1] + v[2]*p[2])*(1-Cos(theta)) + p[2]*Cos(theta) + (-v[1]*p[0] + v[0]*p[1])*Sin(theta);
	
	return result;
}


/*
 * Takes as input a normalized array. Returns a normalized perpendicular array.
 */
array<double, 3> GetPerp(array<double, 3> v)
{
	// Pick m with v[m] != 0
	int m = -1;
        if(v[0] != 0) m = 0;
        else if(v[1] != 0) m = 1;
        else m = 2;

	// Pick n with n != m
	int n = -1;
        if (m == 0) n = 1;
        else n = 0;

	// Create result
	array<double, 3> result = {{0., 0., 0.}};
	result[n] = v[m] / TMath::Sqrt(v[m]*v[m] + v[n]*v[n]);
        result[m] = -v[n] / TMath::Sqrt(v[m]*v[m] + v[n]*v[n]);

	return result;
}


/*
 * Returns the fraction of cone volume is contained within the TPC volume
 */
double TrueConeContainment(double *pos_, double *dir_, double len, double theta, double ppcm = -1)
{
	// Set ppcm to tuned value if nothing is given as input
	if (ppcm == -1)
	{
		int Npoints = 63;
        	double volume = Power(len, 3) * Power(Tan(theta), 2) * Pi()/3;
        	ppcm = (double)Npoints / Power(volume, (1./3.));
	}

	// Reset random number seed based on time
	srand (time(NULL));	
	
	// Create std::array out of input array pointers
	array<double, 3> pos = {{pos_[0], pos_[1], pos_[2]}};
	array<double, 3> dir = {{dir_[0], dir_[1], dir_[2]}};

	// Determine number of segments along length of cone
	int Nlong = (int)len*ppcm + 1;
	
	// Create holder for all the volume points 
	std::vector< array<double, 3> > *volPoints = new std::vector< array<double, 3> >(0);

	// Generate a set of points uniform in volume
	// Iterate down the length of the cone
	for (int i = 0; i < Nlong; i++)
	{
		// How far down the length we are
		double len_i = (double)i/ppcm;
		// Compute number of points in the radial direction
		int Nradial = (int)(len_i * Tan(theta) * ppcm) + 1;

		// Iterate over radii
		for (int j = 0; j < Nradial; j++)
		{
			// On central axis
			if (j == 0)
			{	
				double x = pos[0] + dir[0]*len_i;
				double y = pos[1] + dir[1]*len_i;
				double z = pos[2] + dir[2]*len_i;
		
				array<double, 3> p { {x, y, z } };
				volPoints->push_back(p);
			}

			// Off central axis
			if (j > 0)
			{
				// Moving outward in radius
				double radius_j = len_i * ((double)j/(double)(Nradial-1)) * Tan(theta);
				// Number of points in the circumference of radius radius_j
				double Ncirc_j = (int)(2 * Pi() * radius_j * ppcm);
			
				// Set random angle offset of initial perpendicular vector	
				int offset = 2 * Pi() * (double)(rand()%99) / 100.;
                                array<double, 3> perp = GetPerp(dir);
				perp = RotatePoint(perp, dir, offset);
	
				// Compute dtheta going around the cone central axis, loop over
				double dtheta = 2 * Pi() / (double)Ncirc_j;
				for (int k = 0; k < Ncirc_j; k++)
				{
					array<double, 3> rotPerp = RotatePoint(perp, dir, dtheta*(double)k);
					double x = pos[0] + dir[0]*len_i + rotPerp[0]*radius_j;
					double y = pos[1] + dir[1]*len_i + rotPerp[1]*radius_j;
					double z = pos[2] + dir[2]*len_i + rotPerp[2]*radius_j;
				
					array<double, 3> p { {x, y, z} };
					volPoints->push_back(p);	
				}
			}
		}	
	}

	// Loop through arrays in volPoints and count fraction within fiducial volume
	int Npoints = volPoints->size();
	int NTPC = 0;
	for (array<double, 3> q : *volPoints)
		if (PointInTPCVolume(q)) NTPC ++;

	//cout<<"Total nuber of points: "<<Npoints<<endl;
	//cout<<"Contained points: "<<NTPC<<endl;
	volPoints->clear();
	delete volPoints;

	return (double)NTPC/(double)Npoints;
}


/*
 * Does check on TPC volume containment: create cone that slowly
 * protrudes into the TPC volume face-first
 */
void CheckContainment(double shower_length = 250., double shower_open_angle = .25)
{	
	double pos[] = {125., 0., -shower_length -1.01};
	double dir[] = {0., 0., 1.};
	int Npoints = 63;
	double volume = Power(shower_length, 3) * Power(Tan(shower_open_angle), 2) * Pi()/3;
	cout<<"Volume: "<<volume<<endl;
	cout<<"Length: "<<Power(volume, (1./3.))<<endl;
	double ppcm = (double)Npoints / Power(volume, (1./3.));
	cout<<ppcm<<endl;

	int N = 1000;
	double z[N], fidVol_true[N],  fidVol_app[N];
	double Dmax = 0.;
	double RMS = 0.;

	clock_t t = clock();
	for (int i = 0; i < N; i++)
	{
		pos[2] = pos[2] + (shower_length+1)/(double)N;
		z[i] = pos[2];
		fidVol_true[i] = (Power(shower_length, 3) - Power(-pos[2], 3)) / Power(shower_length, 3);
		fidVol_app[i] = TrueConeContainment(pos, dir, shower_length, shower_open_angle, ppcm); 
		if (Abs(fidVol_true[i]-fidVol_app[i]) > Dmax && pos[2] >= -shower_length) Dmax = Abs(fidVol_true[i]-fidVol_app[i]);
		if (pos[2] >= -shower_length) RMS += Power(fidVol_true[i]-fidVol_app[i],2)/(double)N;
	}
	t = clock() - t;
	RMS = Sqrt(RMS);
	cout<<Form("Run time for ppcm = %.2f: %f milliseconds", ppcm, (double)t/(double)CLOCKS_PER_SEC/(double)N)<<endl;
	
	TGraph *trueVol = new TGraph(N, z, fidVol_true);
	TGraph *appVol = new TGraph(N, z, fidVol_app);
	
	trueVol->SetLineColor(kBlack);
	trueVol->SetLineWidth(3);
	appVol->SetLineColor(kRed);
	appVol->SetLineWidth(3);
	
	TLegend *leg = new TLegend(0.55, 0.15, .85, .3);
	leg->AddEntry(trueVol, "True Fid Vol", "l");
	leg->AddEntry(appVol, "Approx Fid Vol", "l");
	leg->SetBorderSize(0);

	TPaveText *pt = new TPaveText(0.55, .3, .85, .65, "blNDC");
	pt->SetBorderSize(0);
	pt->SetFillColor(0);
	pt->SetTextAlign(12);
	pt->AddText(Form("#Delta_{max} = %.2f%%", Dmax*100.));
	pt->AddText(Form("#Delta_{RMS} = %.2f%%", RMS*100.));
	pt->AddText("");
	pt->AddText(Form("L_{cone} = %.0f cm", shower_length));
	pt->AddText(Form("#theta_{cone} = %.2f", shower_open_angle));
	pt->AddText(Form("%.1f points/cm", ppcm));

	TMultiGraph *mg = new TMultiGraph();
	mg->SetTitle("Fid Vol Approx. Check;Z of Cone Vertex (cm); Fraction Fiducial Volume");
	mg->Add(trueVol);
	mg->Add(appVol);
	
	TCanvas *can = new TCanvas("can", "Fiducial volume approximation check", 200, 10, 1400, 1000);
	mg->Draw("AC");
	can->Modified(); can->Update();
	mg->GetXaxis()->SetRangeUser(-shower_length, 0.);
	mg->GetYaxis()->SetRangeUser(0., 1.);
	leg->Draw();
	pt->Draw();
	can->Print(Form("plots/FidVolApproximationCheck_%.2fcm_%.2fangle.png", shower_length, shower_open_angle));

	delete trueVol;
	delete appVol;
	delete leg;
	delete pt;
	delete mg;
	delete can;
}


/*
 * Make a new TTree in new root file that contains processed data from the analyzer output
 */
void MakeCompressedTTree()
{
	int N_events = 0;
	TChain *MCsample = new TChain("MCsample");
        ifstream ifs("nue.txt");
        string fileName;
        while(std::getline(ifs,fileName))
        {
                fileName = fileName + "/robertoana/pandoratree";
                if(MCsample->Add(fileName.c_str())==0) cout<<"Failed to add file "<<fileName<<endl;
                else N_events++;
        }
        cout<<"TChain "<<MCsample->GetName()<<" has been created"<<endl;

        // Get number of entries and create counter
        N_events *= 500;
        int n_event = 0;
        cout<<"Number of events in chain: "<<N_events<<endl;

        // Create TTreeReader from this chain
        TTreeReader reader(MCsample);

        // Create TTreeReaderValues for the relevant branches
        // True shower energy (fiducial) obtained by looping through all daughters
        TTreeReaderValue<double> true_E(reader, "true_shower_E_dep");
        TTreeReaderValue<double> true_E_fid(reader, "true_shower_E_dep_fiducial");

        // Reconstructed shower start vertex
        TTreeReaderValue< std::vector<double> > start_x(reader, "shower_start_x");
        TTreeReaderValue< std::vector<double> > start_y(reader, "shower_start_y");
        TTreeReaderValue< std::vector<double> > start_z(reader, "shower_start_z");

        // Reconstructed shower start direction
        TTreeReaderValue< std::vector<double> > dir_x(reader, "shower_dir_x");
        TTreeReaderValue< std::vector<double> > dir_y(reader, "shower_dir_y");
        TTreeReaderValue< std::vector<double> > dir_z(reader, "shower_dir_z");

        // PDG code of matched showers and energy of showers in 3 planes
        TTreeReaderValue< std::vector<int> > matched_showers(reader, "matched_showers");
        TTreeReaderValue< std::vector<double> > matched_showers_energy(reader, "matched_showers_energy");
        TTreeReaderValue< std::vector< std::vector<double> > > shower_energy(reader, "shower_energy");
        TTreeReaderValue< std::vector<double> > shower_length(reader, "shower_length");
        TTreeReaderValue< std::vector<double> > shower_open_angle(reader, "shower_open_angle");

        // Energy and PDG of true neutrino daughters
        TTreeReaderValue< std::vector<double> > daughters_E(reader, "nu_daughters_E");
        TTreeReaderValue< std::vector<int> > daughters_pdg(reader, "nu_daughters_pdg");

        // Vertex and momentum of neutrino daughters
        TTreeReaderValue< std::vector<double> > daughters_vx(reader, "nu_daughters_vx");
        TTreeReaderValue< std::vector<double> > daughters_vy(reader, "nu_daughters_vy");
        TTreeReaderValue< std::vector<double> > daughters_vz(reader, "nu_daughters_vz");

        // Vertex and momentum of neutrino daughters
        TTreeReaderValue< std::vector<double> > daughters_px(reader, "nu_daughters_px");
        TTreeReaderValue< std::vector<double> > daughters_py(reader, "nu_daughters_py");
        TTreeReaderValue< std::vector<double> > daughters_pz(reader, "nu_daughters_pz");

        // Event, subrun, run
        TTreeReaderValue<unsigned int> event(reader, "event");
        TTreeReaderValue<unsigned int> subrun(reader, "subrun");
        TTreeReaderValue<unsigned int> run(reader, "run");

	// Make output file and tree
	TFile *outFile = new TFile("CompressedPandoraLEEShower.root", "UPDATE");
	TTree *outTree = new TTree("ShowerAnalysisTree","Compressed Tree for Shower Analysis");
	
	// Make branch references
	double _true_shower_E = -1;
	double _true_shower_E_fiducial = -1;
	double _shower_x = -1;
	double _shower_y = -1;
	double _shower_z = -1;
	double _shower_dir_x = -1;
	double _shower_dir_y = -1;
	double _shower_dir_z = -1;
	int _num_showers = -1;
	double _shower_energy = -1;
	double _shower_length = -1;
	double _shower_open_angle = -1;
	double _true_primary_e_energy = -1;
	double _shower_containment = -1;

	// Make branch addresses
	outTree->Branch("true_shower_E",&_true_shower_E,"true_shower_E/d");
	outTree->Branch("true_shower_E_fiducial",&_true_shower_E_fiducial,"true_shower_E_fiducial/d");
	outTree->Branch("shower_x",&_shower_x,"shower_x/d");
	outTree->Branch("shower_y",&_shower_y,"shower_y/d");
	outTree->Branch("shower_z",&_shower_z,"shower_z/d");
	outTree->Branch("shower_dir_x",&_shower_dir_x,"shower_dir_x/d");
	outTree->Branch("shower_dir_y",&_shower_dir_y,"shower_dir_y/d");
	outTree->Branch("shower_dir_z",&_shower_dir_z,"shower_dir_z/d");
	outTree->Branch("num_showers",&_num_showers,"num_showers/i");
	outTree->Branch("shower_energy",&_shower_energy,"shower_energy/d");
	outTree->Branch("shower_length",&_shower_length,"shower_length/d");
	outTree->Branch("shower_open_angle",&_shower_open_angle,"shower_open_angle/d");
	outTree->Branch("true_primary_e_energy",&_true_primary_e_energy,"true_primary_e_energy/d");
	outTree->Branch("shower_containment",&_shower_containment,"shower_containment/d");

	int Naccepted = 0;

	while(reader.Next())
        {
                // Display progress
                if (n_event % 5000 == 0)
                        cout<<(double)n_event/(double)N_events*100<<"\%"<<endl;
                n_event++;
		
		_true_shower_E = *true_E;
		_true_shower_E_fiducial = *true_E_fid;

                // Compute maximum electron energy and difference with true shower energy       
                int i_e_mcp = -1;
                double E_max = 0;
                for (int i = 0; i < daughters_E->size(); i++)
                        if (TMath::Abs(daughters_pdg->at(i)) == 11 && daughters_E->at(i) > E_max)
                        {
                                E_max = daughters_E->at(i);
                                i_e_mcp = i;
                        }

                // Skip event if no electron is found
                if(i_e_mcp == -1) continue;

                // Find index of primary electron-matched shower with highest energy
                std::vector<int> i_e_vec = {};
                for (int i = 0; i < matched_showers->size(); i++)
                {
                        if (matched_showers->at(i) == 11 && Abs(matched_showers_energy->at(i) - daughters_E->at(i_e_mcp)) < 0.001)
                                i_e_vec.push_back(i);
                }

                // Skip event if no shower matched to electron
                if (i_e_vec.size() == 0)
                        continue;

		Naccepted++;
		//if(Naccepted > 20) break;

                // Find the shower portion closest to the true electron origin if there are multiple matched showers
                int i_e = i_e_vec.at(0);
                if (i_e_vec.size() > 1)
                {
                        double d_0_x = start_x->at(i_e) - daughters_vx->at(i_e_mcp);
                        double d_0_y = start_y->at(i_e) - daughters_vy->at(i_e_mcp);
                        double d_0_z = start_z->at(i_e) - daughters_vz->at(i_e_mcp);
                        double d_0 = Sqrt( d_0_x*d_0_x + d_0_y*d_0_y + d_0_z*d_0_z);
                        for (int i = 1; i < i_e_vec.size(); i++)
                        {
                                d_0_x = start_x->at(i_e_vec.at(i)) - daughters_vx->at(i_e_mcp);
                                d_0_y = start_y->at(i_e_vec.at(i)) - daughters_vy->at(i_e_mcp);
                                d_0_z = start_z->at(i_e_vec.at(i)) - daughters_vz->at(i_e_mcp);
                                if (Sqrt( d_0_x*d_0_x + d_0_y*d_0_y + d_0_z*d_0_z) < d_0)
                                {
                                        d_0 = Sqrt( d_0_x*d_0_x + d_0_y*d_0_y + d_0_z*d_0_z);
                                        i_e = i_e_vec.at(i);
                                }
                        }
                }

		_shower_x = start_x->at(i_e);
		_shower_y = start_y->at(i_e);
		_shower_z = start_z->at(i_e);

		_num_showers = i_e_vec.size();
						
		// Reco shower direction
                double x_dir, y_dir, z_dir;
                double dir_norm;
                double shower_half_angle;

		x_dir = dir_x->at(i_e);
                y_dir = dir_y->at(i_e);
                z_dir = dir_z->at(i_e);

                dir_norm = Power(x_dir,2) + Power(y_dir,2) + Power(z_dir,2);
                dir_norm = TMath::Sqrt(dir_norm);
	
		_shower_dir_x = x_dir/dir_norm;
                _shower_dir_y = y_dir/dir_norm;
                _shower_dir_z = z_dir/dir_norm;

                //if (Abs(E_max-*true_E) > .01) continue;

                std::vector<double> percRecoShowerFid = {};
                std::vector<double> volume = {};
                for (int i = 0; i < i_e_vec.size(); i++)
                {
                        x_dir = dir_x->at(i_e_vec.at(i));
                        y_dir = dir_y->at(i_e_vec.at(i));
                        z_dir = dir_z->at(i_e_vec.at(i));

                        dir_norm = Power(x_dir,2) + Power(y_dir,2) + Power(z_dir,2);
                        dir_norm = TMath::Sqrt(dir_norm);

                        double shower_dir[3] = {x_dir/dir_norm, y_dir/dir_norm, z_dir/dir_norm};
                        double shower_start[3] = {start_x->at(i_e_vec.at(i)), start_y->at(i_e_vec.at(i)), start_z->at(i_e_vec.at(i))};
                        shower_half_angle = shower_open_angle->at(i_e_vec.at(i)) / 2;

                        percRecoShowerFid.push_back(TrueConeContainment(shower_start, shower_dir, shower_length->at(i_e_vec.at(i)), shower_half_angle));
                        volume.push_back(Pi()/3 * Power(shower_length->at(i_e_vec.at(i)), 3) * Power(Tan(shower_half_angle),2));
                }

                double volContained = 0;
                double vol = 0;
                for ( int i = 0; i < i_e_vec.size(); i++)
                {
                        volContained += percRecoShowerFid.at(i) * volume.at(i);
                        vol += volume.at(i);
                }
                double totalShowerContainment = volContained / vol;

                double shower_E = 0;
                for(int i = 0; i < i_e_vec.size(); i++) shower_E += shower_energy->at(i_e_vec.at(i)).at(2);
		
		_shower_energy = shower_E;
		_shower_length = shower_length->at(i_e);
		_shower_open_angle = shower_open_angle->at(i_e);
		_true_primary_e_energy = E_max;
		_shower_containment = totalShowerContainment;
		
		outTree->Fill();
	}

	outTree->Write();
	outFile->Close();
	delete outFile;
	delete MCsample;
}


/*
 * Create plots of variables and compare between truth level containment/non-containment
 */
void MakeContainmentComparisonPlots(string output_prefix = "")
{
	// If output_prefix is not empty, append underscore after
	if (output_prefix.length() > 0)
		output_prefix.append("_");

	// Set global variables for ROOT
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetLegendBorderSize(0);
	
	TFile *f = new TFile("CompressedPandoraLEEShower.root", "READ");
	TTree* MCsample = (TTree*)f->Get("ShowerAnalysisTree");

	// Create chain of ROOT file samples by reading in all files in nue.txt
	Long64_t N_events = MCsample->GetEntries();
        cout<<MCsample->GetName()<<" has "<<N_events<<" events"<<endl;	

	// Create TTreeReader from this chain
	TTreeReader reader(MCsample);

	// Create TTreeReaderValues for the relevant branches
	// True shower energy (fiducial) obtained by looping through all daughters
	TTreeReaderValue<double> true_E(reader, "true_shower_E");
	TTreeReaderValue<double> true_E_fid(reader, "true_shower_E_fiducial");
	
	// Reconstructed shower start vertex
	TTreeReaderValue<double> start_x(reader, "shower_x");
	TTreeReaderValue<double> start_y(reader, "shower_y");
	TTreeReaderValue<double> start_z(reader, "shower_z");	
		
	// Reconstructed shower start direction
	TTreeReaderValue<double> dir_x(reader, "shower_dir_x");
	TTreeReaderValue<double> dir_y(reader, "shower_dir_y");
	TTreeReaderValue<double> dir_z(reader, "shower_dir_z");

	// Number of matched showers, total energy in Y plane length and open angle
	// of shower closest to true vertex
	TTreeReaderValue<unsigned int> num_showers(reader, "num_showers");
	TTreeReaderValue<double> shower_energy(reader, "shower_energy");
	TTreeReaderValue<double> shower_length(reader, "shower_length");
	TTreeReaderValue<double> shower_open_angle(reader, "shower_open_angle");

	// True primary electron energy and shower containment
	TTreeReaderValue<double> true_primary_E(reader, "true_primary_e_energy");
	TTreeReaderValue<double> shower_containment(reader, "shower_containment");
	
	// Number of showers matched to electrons
	TH1F *TH1_numShowers = new TH1F("TH1_numShowers", "", 10, 0., 5);
	TH1_numShowers->SetTitle("; Number ofEnergy-Electron-Matched Showers;");
 
	// Make histogram to store energy deposition vs true electron daughter energy
	TH2F *TH2_trueE = new TH2F("TH2_trueE", "", 50, 0, 5, 50, 0, 5);
	TH2_trueE->SetTitle(";True Max E_{e} Daughter (GeV);True Total Shower Energy (GeV)");
	
	// Make histogram to store true energy deposition vs reco shower energy
	TH2F *TH2_recoTrueE = new TH2F("TH2_recoTrueE", "", 50, 0, 5, 50, 0, 5);
	TH2_recoTrueE->SetTitle(";Total Reco Shower Energy (GeV); True Shower Energy (GeV)");

	// Make histogram for containment
	TH2F* TH2_trueRecoFid = new TH2F("TH2_trueRecoFid", "", 50, 0, 1, 50, 0, 1);
	TH2_trueRecoFid->SetTitle(";Fiducial Fraction of True Shower Energy; Fiducial Fraction of Reco Shower Volume");

	// Histogram of absolute energy difference
	TH1F* TH1_checkTrueEDiff= new TH1F("TH1_checkTrueEDiff", "", 50, 0, 50);
	TH1_checkTrueEDiff->SetTitle(";|E_{shower}-E_{e^{-}}| (MeV);");
	TH1_checkTrueEDiff->SetLineWidth(2);

	// Histogram of primary electron vertices vs true containment
	TH2F* TH2_vxTrueFid = new TH2F("TH2_vxTrueFid", "", 50, .8, 1, 50, 0, 256.35);
	TH2_vxTrueFid->SetTitle(";True Shower Containment (%);Primary Electron X-Vertex (cm)");

	TH2F* TH2_vyTrueFid = new TH2F("TH2_vyTrueFid", "", 50, .8, 1, 50, -116.5, 116.5);
        TH2_vyTrueFid->SetTitle(";True Shower Containment (%);Primary Electron Y-Vertex (cm)");

	TH2F* TH2_vzTrueFid = new TH2F("TH2_vzTrueFid", "", 50, .8, 1, 50, 0, 1036.8);
        TH2_vzTrueFid->SetTitle(";True Shower Containment (%);Primary Electron Z-Vertex (cm)");
	
	// 1D Histograms of Reco Shower start positions
	TH1* TH1_ShowerZpass = new TH1F("TH1_ShowerZpass", "", 50, 0, 1036.8);
	TH1_ShowerZpass->SetTitle(";Reco Shower Start-Z (cm);");
	TH1_ShowerZpass->SetLineWidth(2);
	
	TH1* TH1_ShowerZfail = new TH1F("TH1_ShowerZfail", "", 50, 0, 1036.8);
        TH1_ShowerZfail->SetTitle(";Reco Shower Start-Z (cm);");
	TH1_ShowerZfail->SetLineColor(kRed+2);
	TH1_ShowerZfail->SetLineWidth(2);
	
	TH1* TH1_ShowerYpass = new TH1F("TH1_ShowerYpass", "", 50, -116.5, 116.5);
        TH1_ShowerYpass->SetTitle(";Reco Shower Start-Y (cm);");
        TH1_ShowerYpass->SetLineWidth(2);

        TH1* TH1_ShowerYfail = new TH1F("TH1_ShowerYfail", "", 50, -116.5, 116.5);
        TH1_ShowerYfail->SetTitle(";Reco Shower Start-Y (cm);");
        TH1_ShowerYfail->SetLineColor(kRed+2);
        TH1_ShowerYfail->SetLineWidth(2);

	TH1* TH1_ShowerXpass = new TH1F("TH1_ShowerXpass", "", 50, 0, 256.35);
        TH1_ShowerXpass->SetTitle(";Reco Shower Start-X (cm);");
        TH1_ShowerXpass->SetLineWidth(2);

        TH1* TH1_ShowerXfail = new TH1F("TH1_ShowerXfail", "", 50, 0, 256.35);
        TH1_ShowerXfail->SetTitle(";Reco Shower Start-X (cm);");
        TH1_ShowerXfail->SetLineColor(kRed+2);
        TH1_ShowerXfail->SetLineWidth(2);

	// 1D Histograms for Reco Shower open angle
	TH1* TH1_ShowerOpenAnglePass = new TH1F("TH1_ShowerOpenAnglePass", "", 50, 0, Pi()/2);
        TH1_ShowerOpenAnglePass->SetTitle(";Reco Shower Open Angle;");
        TH1_ShowerOpenAnglePass->SetLineWidth(2);

        TH1* TH1_ShowerOpenAngleFail = new TH1F("TH1_ShowerOpenAngleFail", "", 50, 0, Pi()/2);
        TH1_ShowerOpenAngleFail->SetTitle(";Reco Shower Open Angle;");
        TH1_ShowerOpenAngleFail->SetLineColor(kRed+2);
        TH1_ShowerOpenAngleFail->SetLineWidth(2);

	// 1D Histograms for Reco Shower energy
        TH1* TH1_ShowerEnergyPass = new TH1F("TH1_ShowerEnergyPass", "", 50, 0, 4);
        TH1_ShowerEnergyPass->SetTitle(";Reco Shower Energy (GeV);");
        TH1_ShowerEnergyPass->SetLineWidth(2);

        TH1* TH1_ShowerEnergyFail = new TH1F("TH1_ShowerEnergyFail", "", 50, 0, 4);
        TH1_ShowerEnergyFail->SetTitle(";Reco Shower Energy (GeV);");
        TH1_ShowerEnergyFail->SetLineColor(kRed+2);
        TH1_ShowerEnergyFail->SetLineWidth(2);

	// 1D Histograms for Reco Shower length
        TH1* TH1_ShowerLengthPass = new TH1F("TH1_ShowerLengthPass", "", 50, 0, 450);
        TH1_ShowerLengthPass->SetTitle(";Reco Shower Length (cm);");
        TH1_ShowerLengthPass->SetLineWidth(2);

        TH1* TH1_ShowerLengthFail = new TH1F("TH1_ShowerLengthFail", "", 50, 0, 450);
        TH1_ShowerLengthFail->SetTitle(";Reco Shower Length (cm);");
        TH1_ShowerLengthFail->SetLineColor(kRed+2);
        TH1_ShowerLengthFail->SetLineWidth(2);

	// 1D Histogram for Reco Containment
	TH1* TH1_ShowerContainmentPass = new TH1F("TH1_ShowerContainmentPass", "", 50, 0, 1);
        TH1_ShowerContainmentPass->SetTitle(";Reco Shower Containment;");
        TH1_ShowerContainmentPass->SetLineWidth(2);

        TH1* TH1_ShowerContainmentFail = new TH1F("TH1_ShowerContainmentFail", "", 50, 0, 1);
        TH1_ShowerContainmentFail->SetTitle(";Reco Shower Containment;");
        TH1_ShowerContainmentFail->SetLineColor(kRed+2);
        TH1_ShowerContainmentFail->SetLineWidth(2);

	// 2D Histogram for Positions
	// XY
	TH2* TH2_ShowerXYpass = new TH2F("TH2_ShowerXYpass", "", 50, 0, 256.35, 50, -116.5, 116.5);
	TH2_ShowerXYpass->SetTitle("Pass;Reco Shower Start-X (cm);Reco Shower Start-Y (cm)");
	TH2_ShowerXYpass->SetMarkerStyle(20);
	TH2_ShowerXYpass->SetMarkerSize(0.35);
	TH2_ShowerXYpass->SetMarkerColorAlpha(kGreen, .4);

	TH2* TH2_ShowerXYfail = new TH2F("TH2_ShowerXYfail", "", 50, 0, 256.35, 50, -116.5, 116.5);
        TH2_ShowerXYfail->SetTitle("Fail;Reco Shower Start-X (cm);Reco Shower Start-Y (cm)");
	TH2_ShowerXYfail->SetMarkerStyle(20);
	TH2_ShowerXYfail->SetMarkerSize(0.35);
	TH2_ShowerXYfail->SetMarkerColorAlpha(kRed, .4);

	// XZ
	TH2* TH2_ShowerXZpass = new TH2F("TH2_ShowerXZpass", "", 50, 0, 256.35, 50, 0, 1036.8);
        TH2_ShowerXZpass->SetTitle("Pass;Reco Shower Start-X (cm);Reco Shower Start-Z (cm)");
        TH2_ShowerXZpass->SetMarkerStyle(20);
        TH2_ShowerXZpass->SetMarkerSize(0.35);
        TH2_ShowerXZpass->SetMarkerColorAlpha(kGreen, .4);

        TH2* TH2_ShowerXZfail = new TH2F("TH2_ShowerXZfail", "", 50, 0, 256.35, 50, 0, 1036.8);
        TH2_ShowerXZfail->SetTitle("Fail;Reco Shower Start-X (cm);Reco Shower Start-Z (cm)");
        TH2_ShowerXZfail->SetMarkerStyle(20);
        TH2_ShowerXZfail->SetMarkerSize(0.35);
        TH2_ShowerXZfail->SetMarkerColorAlpha(kRed, .4);

	// YZ
	TH2* TH2_ShowerYZpass = new TH2F("TH2_ShowerYZpass", "", 50, -116.5, 116.5, 50, 0, 1036.8);
        TH2_ShowerYZpass->SetTitle("Pass;Reco Shower Start-Y (cm);Reco Shower Start-Z (cm)");
        TH2_ShowerYZpass->SetMarkerStyle(20);
        TH2_ShowerYZpass->SetMarkerSize(0.35);
        TH2_ShowerYZpass->SetMarkerColorAlpha(kGreen, .4);

        TH2* TH2_ShowerYZfail = new TH2F("TH2_ShowerYZfail", "", 50, -116.5, 116.5, 50, 0, 1036.8);
        TH2_ShowerYZfail->SetTitle("Fail;Reco Shower Start-Y (cm);Reco Shower Start-Z (cm)");
        TH2_ShowerYZfail->SetMarkerStyle(20);
        TH2_ShowerYZfail->SetMarkerSize(0.35);
        TH2_ShowerYZfail->SetMarkerColorAlpha(kRed, .4);
	
	/*
	TH1* TH1_TrueShowerContainmentPass = new TH1F("TH1_TrueShowerContainmentPass","",50,0,1);
	TH1_TrueShowerContainmentPass->SetTitle(";True Shower Containment Fraction;");
	
	TH1* TH1_TrueShowerContainmentFail = new TH1F("TH1_TrueShowerContainmentFail","",50,0,1);
        TH1_TrueShowerContainmentFail->SetTitle(";True Shower Containment Fraction;");
	TH1_TrueShowerContainmentFail->SetLineColor(kRed+2);	
	*/

	// Count total number of events accepted
	double NeventsAccepted = 0;

	while(reader.Next())
	{
		TH1_numShowers->Fill(*num_showers);
	
		TH1_checkTrueEDiff->Fill(Abs(*true_primary_E-*true_E)*1000.);	
	
		// Skip if true shower differs from true electron energy by more than 10 MeV
		if (Abs(*true_primary_E-*true_E) > .01) continue;
		
		NeventsAccepted++;
		
		// Fill 2D histogram to compare true electron energy to true shower energy
		TH2_trueE->Fill(*true_primary_E, *true_E );

		TH2_recoTrueE->Fill(*shower_energy, *true_E);
		TH2_trueRecoFid->Fill( (*true_E_fid)/(*true_E), *shower_containment);
	
		// Desired minimum true shower containment	
		double cut_true_shower_containment = 0.8;
	
		double cut_containment_agreement = 0.05;
		double trueContainment = (*true_E_fid)/(*true_E);
		if (Abs(trueContainment-*shower_containment)/trueContainment <= cut_containment_agreement)
		//if (trueContainment >= cut_true_shower_containment)
		{
			TH1_ShowerXpass->Fill(*start_x);
			TH1_ShowerYpass->Fill(*start_y);
			TH1_ShowerZpass->Fill(*start_z);
			TH1_ShowerOpenAnglePass->Fill(*shower_open_angle);
			TH1_ShowerContainmentPass->Fill(*shower_containment);
			TH1_ShowerLengthPass->Fill(*shower_length);
			TH1_ShowerEnergyPass->Fill(*shower_energy);
	
			TH2_ShowerXYpass->Fill(*start_x, *start_y);
			TH2_ShowerXZpass->Fill(*start_x, *start_z);
			TH2_ShowerYZpass->Fill(*start_y, *start_z);
		}
	
		else
		{
			TH1_ShowerXfail->Fill(*start_x);
			TH1_ShowerYfail->Fill(*start_y);
			TH1_ShowerZfail->Fill(*start_z);
			TH1_ShowerOpenAngleFail->Fill(*shower_open_angle);
			TH1_ShowerContainmentFail->Fill(*shower_containment);
			TH1_ShowerLengthFail->Fill(*shower_length);
			TH1_ShowerEnergyFail->Fill(*shower_energy);

			TH2_ShowerXYfail->Fill(*start_x, *start_y);
			TH2_ShowerXZfail->Fill(*start_x, *start_z);
			TH2_ShowerYZfail->Fill(*start_y, *start_z);
		}
	}

	TCanvas *can = new TCanvas("can", "can", 800, 800);

	TH2_trueRecoFid->Draw("COLZ");
	can->Print( Form("plots/%sTH2_TrueShowerContainment_vs_RecoShowerContainment.png", output_prefix.c_str()) );

	TH1_ShowerXfail->Draw("hist");
	TH1_ShowerXpass->Draw("hist same");
	can->Print( Form("plots/%sTH1_ShowerX.png", output_prefix.c_str()) );
	
	TH1_ShowerYfail->Draw("hist");
        TH1_ShowerYpass->Draw("hist same");
        can->Print( Form("plots/%sTH1_ShowerY.png", output_prefix.c_str()) );

	TH1_ShowerZpass->Draw("hist");
        TH1_ShowerZfail->Draw("hist same");
        can->Print( Form("plots/%sTH1_ShowerZ.png", output_prefix.c_str()) );

	TH1_ShowerEnergyFail->Draw("hist");
        TH1_ShowerEnergyPass->Draw("hist same");
        can->Print( Form("plots/%sTH1_ShowerEnergy.png", output_prefix.c_str()) );

	TH1_ShowerLengthFail->Draw("hist");
        TH1_ShowerLengthPass->Draw("hist same");
        can->Print( Form("plots/%sTH1_ShowerLength.png", output_prefix.c_str()) );

	TH1_ShowerOpenAnglePass->Draw("hist");
        TH1_ShowerOpenAngleFail->Draw("hist same");
        can->Print( Form("plots/%sTH1_ShowerOpenAngle.png", output_prefix.c_str()) );

	TH1_ShowerContainmentPass->Draw("hist");
        TH1_ShowerContainmentFail->Draw("hist same");
        can->Print( Form("plots/%sTH1_ShowerContainment.png", output_prefix.c_str()) );

	TH2_ShowerXYpass->Draw("SCAT");
	TH2_ShowerXYfail->Draw("SCAT same");
	can->Print( Form("plots/%sTH2_ShowerXY_scatter.png", output_prefix.c_str()) );
	
	TH2_ShowerXZpass->Draw("SCAT");
        TH2_ShowerXZfail->Draw("SCAT same");
        can->Print( Form("plots/%sTH2_ShowerXZ_scatter.png", output_prefix.c_str()) );

	TH2_ShowerYZpass->Draw("SCAT");
        TH2_ShowerYZfail->Draw("SCAT same");
        can->Print( Form("plots/%sTH2_ShowerYZ_scatter.png", output_prefix.c_str()) );


	TCanvas *can2 = new TCanvas("can2","can",1600,800);
        can2->Divide(2,1);

	can2->cd(1);
	TH2_ShowerXYpass->Draw("COLZ");
	can2->cd(2);
	TH2_ShowerXYfail->Draw("COLZ");
	can2->Print( Form("plots/%sTH2_ShowerXY_colz.png", output_prefix.c_str()) );

	can2->cd(1);
        TH2_ShowerXZpass->Draw("COLZ");
        can2->cd(2);
        TH2_ShowerXZfail->Draw("COLZ");
        can2->Print( Form("plots/%sTH2_ShowerXZ_colz.png", output_prefix.c_str()) );

	can2->cd(1);
        TH2_ShowerYZpass->Draw("COLZ");
        can2->cd(2);
        TH2_ShowerYZfail->Draw("COLZ");
        can2->Print( Form("plots/%sTH2_ShowerYZ_colz.png", output_prefix.c_str()) );

	can->cd();
	TH2_trueE->Draw("COLZ");
        can->Print( Form("plots/%sTH2_TrueShowerEnergy_vs_MCParticleEnergy.png", output_prefix.c_str()) );	
	
	TH2_recoTrueE->Draw("COLZ");
        can->Print( Form("plots/%sTH2_TrueShowerEnergy_vs_RecoShowerEnergy.png", output_prefix.c_str()) );
	
	TH1_checkTrueEDiff->Draw("hist");
	can->Print( Form("plots/%sTH1_AbsoluteDifferenceOfShowerAndElectronEnergy.png", output_prefix.c_str()) );
	
	TH2_trueRecoFid->Draw("COLZ");
	can->Print( Form("plots/%sTH2_TrueShowerContainment_vs_RecoShowerContainment.png", output_prefix.c_str()) );	
	
	TH1_numShowers->Draw("hist");
	can->Print( Form("plots/%sTH1_NumberOfMatchedShowers.png", output_prefix.c_str()) );

	delete MCsample;
	delete TH2_trueE;
	delete TH2_recoTrueE;
	delete TH2_trueRecoFid;
	delete TH1_numShowers;
	delete TH1_checkTrueEDiff;
	delete TH2_vxTrueFid;
	delete TH2_vyTrueFid;
	delete TH2_vzTrueFid;
	delete TH1_ShowerXpass;
	delete TH1_ShowerYpass;
	delete TH1_ShowerZpass;
	delete TH1_ShowerOpenAnglePass;
	delete TH1_ShowerContainmentPass;
	delete TH1_ShowerEnergyPass;
	delete TH1_ShowerLengthPass;
	delete TH1_ShowerXfail;
	delete TH1_ShowerYfail;
	delete TH1_ShowerZfail;
	delete TH1_ShowerOpenAngleFail;
	delete TH1_ShowerLengthFail;
	delete TH1_ShowerEnergyFail;
	delete TH1_ShowerContainmentFail;
	delete can;
}


/*
 * Do an ensemble of selections with variable threshold for true shower containment
 */
void TestShowerSelectionEnsemble(string output_prefix = "", double true_shower_containment_requirement = 0.95)
{
	// If output_prefix is not empty, append underscore after
	if (output_prefix.length() > 0)
		output_prefix.append("_");

	// Set global variables for ROOT
	gROOT->Reset();
	gStyle->SetOptStat(0);
	gStyle->SetLegendBorderSize(0);
	
	TFile *f = new TFile("CompressedPandoraLEEShower.root", "READ");
	TTree* MCsample = (TTree*)f->Get("ShowerAnalysisTree");

	// Create chain of ROOT file samples by reading in all files in nue.txt
	Long64_t N_events = MCsample->GetEntries();
        cout<<MCsample->GetName()<<" has "<<N_events<<" events"<<endl;	

	// Create TTreeReader from this chain
	TTreeReader reader(MCsample);

	// Create TTreeReaderValues for the relevant branches
	// True shower energy (fiducial) obtained by looping through all daughters
	TTreeReaderValue<double> true_E(reader, "true_shower_E");
	TTreeReaderValue<double> true_E_fid(reader, "true_shower_E_fiducial");
	
	// Reconstructed shower start vertex
	TTreeReaderValue<double> start_x(reader, "shower_x");
	TTreeReaderValue<double> start_y(reader, "shower_y");
	TTreeReaderValue<double> start_z(reader, "shower_z");	
		
	// Reconstructed shower start direction
	TTreeReaderValue<double> dir_x(reader, "shower_dir_x");
	TTreeReaderValue<double> dir_y(reader, "shower_dir_y");
	TTreeReaderValue<double> dir_z(reader, "shower_dir_z");

	// Number of matched showers, total energy in Y plane length and open angle
	// of shower closest to true vertex
	TTreeReaderValue<unsigned int> num_showers(reader, "num_showers");
	TTreeReaderValue<double> shower_energy(reader, "shower_energy");
	TTreeReaderValue<double> shower_length(reader, "shower_length");
	TTreeReaderValue<double> shower_open_angle(reader, "shower_open_angle");

	// True primary electron energy and shower containment
	TTreeReaderValue<double> true_primary_E(reader, "true_primary_e_energy");
	TTreeReaderValue<double> shower_containment(reader, "shower_containment");

	
	// Make output file and tree	
	TFile *outFile = new TFile(Form("%sSelectionOptimizationTree.root", output_prefix.c_str()),"UPDATE");
	TTree *outTree = new TTree("OptimizationTree","Tree Containing Selection Parameters with Results");

	// Create branches of optimizaton tree
	double radius_x = 0;
	double radius_y = 0;
	double z_1 = 0;
	double z_2 = 0;
	double min_shower_containment = 0;
	double min_shower_energy = 0;
        double min_shower_length = 0;
        double max_shower_length = 0;
        double max_shower_open_angle = 0;

	double mean_true_containment = 0;
	double NacceptedEvents = 0;
        double truePositive = 0;
        double falsePositive = 0;
        double trueNegative = 0;
        double falseNegative = 0;
		
	// Make branch addresses
	outTree->Branch("radius_x",&radius_x,"radius_x/d");
	outTree->Branch("radius_y",&radius_y,"radius_y/d");
	outTree->Branch("z_1",&z_1,"z_1/d");
	outTree->Branch("z_2",&z_2,"z_2/d");
	outTree->Branch("min_shower_containment",&min_shower_containment,"min_shower_containment/d");
	//outTree->Branch("min_shower_energy",&min_shower_energy,"min_shower_energy/d");
        //outTree->Branch("min_shower_length",&min_shower_length,"min_shower_length/d");
        //outTree->Branch("max_shower_length",&max_shower_length,"max_shower_length/d");
        //outTree->Branch("max_shower_open_angle",&max_shower_open_angle,"max_shower_open_angle/d");

	outTree->Branch("mean_true_containment",&mean_true_containment,"mean_true_containment/d");
	outTree->Branch("NacceptedEvents",&NacceptedEvents,"NacceptedEvents/d");
	outTree->Branch("truePositive",&truePositive,"truePositive/d");
	outTree->Branch("falsePositive",&falsePositive,"falsePositive/d");
	outTree->Branch("trueNegative",&trueNegative,"trueNegative/d");
	outTree->Branch("falseNegative",&falseNegative,"falseNegative/d");
	
	// Set bounds over which to vary the selection parameters	
	double vradius_x[2] = {95., 100.};
        double vradius_y[2] = {80., 85.};
        double vz_1[2] = {25., 60.};
        double vz_2[2] = {900., 950.};
        double vmin_shower_containment[2] = {0.7, .85};
        double vmin_shower_energy[2] = {0.1, 2.};
        double vmin_shower_length[2] = {60., 125.};
        double vmax_shower_length[2] = {175, 400.};
        double vmax_shower_open_angle[2] = {0.2, 0.8};

	/*double vradius_x[2] = {5., 10.};
        double vradius_y[2] = {5., 25.};
        double vz_1[2] = {25., 200.};
        double vz_2[2] = {800., 870.};
        double vmin_shower_containment[2] = {0.9, 1.};
        double vmin_shower_energy[2] = {0.1475, 2.};
        double vmin_shower_length[2] = {60., 125.};
        double vmax_shower_length[2] = {175, 225.};
        double vmax_shower_open_angle[2] = {0.2, 0.5};
	*/
	
	int num_tests = 100000;
	for (int test = 0; test < num_tests; test++)
	{
		if (test % (int)((double)num_tests/10.) == 0)
		{
			srand (time(NULL));
			cout<<(double)test/(double)num_tests*100<<"\%"<<endl;
		}

		radius_x = vradius_x[0] + (double)(rand()%999)/1000.*(vradius_x[1]-vradius_x[0]);
	     	radius_y = vradius_y[0] + (double)(rand()%999)/1000.*(vradius_y[1]-vradius_y[0]);
        	z_1 = vz_1[0] + (double)(rand()%999)/1000.*(vz_1[1]-vz_1[0]);
       		z_2 = vz_2[0] + (double)(rand()%999)/1000.*(vz_2[1]-vz_2[0]);
        	min_shower_containment = vmin_shower_containment[0] + (double)(rand()%999)/1000.*(vmin_shower_containment[1]-vmin_shower_containment[0]);
		min_shower_energy = vmin_shower_energy[0] + (double)(rand()%999)/1000.*(vmin_shower_energy[1]-vmin_shower_energy[0]);
                min_shower_length = vmin_shower_length[0] + (double)(rand()%999)/1000.*(vmin_shower_length[1]-vmin_shower_length[0]);
                max_shower_length = vmax_shower_length[0] + (double)(rand()%999)/1000.*(vmax_shower_length[1]-vmax_shower_length[0]);
                max_shower_open_angle = vmax_shower_open_angle[0] + (double)(rand()%999)/1000.*(vmax_shower_open_angle[1]-vmax_shower_open_angle[0]);

		mean_true_containment = 0;
        	NacceptedEvents = 0;
        	truePositive = 0;
        	falsePositive = 0;
        	trueNegative = 0;
        	falseNegative = 0;

		while(reader.Next())
		{
			// Skip if true shower differs from true electron energy by more than 10 MeV
			if (Abs(*true_primary_E-*true_E) > .01) continue;
			// Skip events when the shower begins in the dead region of the Y-plane
			if (*start_z >= 675. && *start_z <=775. /*&& true_shower_containment_requirement > .9*/) continue;	

			NacceptedEvents++;
			
			double Q = Power((*start_x-125.)/radius_x,2) + Power(*start_y/radius_y,2);
			if (*start_z >= z_1 &&
                            *start_z <= z_2 &&
                            Q <= 1. &&
                            *shower_containment >= min_shower_containment )
                     	    //*shower_energy >= min_shower_energy && 
                     	    //*shower_length >= min_shower_length && 
                     	    //*shower_length <= max_shower_length && 
                     	    //*shower_open_angle <= max_shower_open_angle)
			{
				if ((*true_E_fid)/(*true_E) >= true_shower_containment_requirement) truePositive++;
				else falsePositive++;
				mean_true_containment += (*true_E_fid)/(*true_E);
			}
			else 
			{
                        	if ((*true_E_fid)/(*true_E) < true_shower_containment_requirement) trueNegative++;
                        	else falseNegative++;
			}	
		}	
	
		mean_true_containment /= (truePositive+falsePositive);

		// Fill tree and restart ttreeReader
		outTree->Fill();
		reader.SetEntry(0);
	}
	outTree->Write();
	outFile->Close();

	delete outFile;
	delete MCsample;
}


/*
 * Make plots for results of ensemble of cuts
 */
void AnalyzeShowerSelectionTree(string output_prefix = "")
{
        // If output_prefix is not empty, append underscore after
        if (output_prefix.length() > 0)
                output_prefix.append("_");

        // Set global variables for ROOT
        gROOT->Reset();
        gStyle->SetOptStat(0);
        gStyle->SetLegendBorderSize(0);

        TFile *f = new TFile("95p_SelectionOptimizationTree.root", "READ");
        TTree* SelectionTree = (TTree*)f->Get("OptimizationTree");

	TH2* TH2_precisionRecall = new TH2F("TH2_precisionRecall","",50,.89,1.,50,0.,.3);
	TH2_precisionRecall->SetTitle(";Precision [TP/(TP+FP)];Recall [TP/(TP+FN)]");
	TH2_precisionRecall->SetMarkerStyle(20);
        TH2_precisionRecall->SetMarkerSize(0.35);
        TH2_precisionRecall->SetMarkerColorAlpha(kBlack, .4);
	TH2_precisionRecall->GetYaxis()->SetTitleOffset(1.4);
	TH2_precisionRecall->GetXaxis()->SetNdivisions(505);

	TH2* TH2_precisionMinShowerContainment = new TH2F("TH2_precisionMinShowerContainment", "", 50,.9,1, 50,.945,.967);
        TH2_precisionMinShowerContainment->SetTitle(";Minimum Shower Containment;Precision");
        TH2_precisionMinShowerContainment->SetMarkerStyle(20);
        TH2_precisionMinShowerContainment->SetMarkerSize(0.5);
        TH2_precisionMinShowerContainment->SetMarkerColorAlpha(kBlack, .35);
        TH2_precisionMinShowerContainment->GetYaxis()->SetTitleOffset(1.6);
        TH2_precisionMinShowerContainment->GetXaxis()->SetNdivisions(505);	

	TH2* TH2_recallMinShowerContainment = new TH2F("TH2_recallMinShowerContainment", "", 50,.9,1, 50,.075,.25);
        TH2_recallMinShowerContainment->SetTitle(";Minimum Shower Containment;Recall");
        TH2_recallMinShowerContainment->SetMarkerStyle(20);
        TH2_recallMinShowerContainment->SetMarkerSize(0.5);
        TH2_recallMinShowerContainment->SetMarkerColorAlpha(kRed, .35);
        TH2_recallMinShowerContainment->GetYaxis()->SetTitleOffset(1.4);
        TH2_recallMinShowerContainment->GetXaxis()->SetNdivisions(505);
	TH2_recallMinShowerContainment->GetYaxis()->SetAxisColor(kRed);
	TH2_recallMinShowerContainment->GetYaxis()->SetLabelColor(kRed);
	TH2_recallMinShowerContainment->GetYaxis()->SetTitleColor(kRed);

	TH2* TH2_radius_xRadius_y = new TH2F("TH2_radius_xRadius_y", "", 50,50,100, 50,35,85);
        TH2_radius_xRadius_y->SetTitle(";Radius X (cm);Radius Y (cm)");
        TH2_radius_xRadius_y->SetMarkerStyle(20);
        TH2_radius_xRadius_y->SetMarkerSize(0.35);
        TH2_radius_xRadius_y->SetMarkerColorAlpha(kGreen, .4);
        TH2_radius_xRadius_y->GetYaxis()->SetTitleOffset(1.4);
        TH2_radius_xRadius_y->GetXaxis()->SetNdivisions(505);

	TH2* TH2_radius_xRadius_yBad = new TH2F("TH2_radius_xRadius_yBad", "", 50,50,100, 50,35,85);
        TH2_radius_xRadius_yBad->SetTitle(";Radius X (cm);Radius Y (cm)");
        TH2_radius_xRadius_yBad->SetMarkerStyle(20);
        TH2_radius_xRadius_yBad->SetMarkerSize(0.35);
        TH2_radius_xRadius_yBad->SetMarkerColorAlpha(kRed, .4);
        TH2_radius_xRadius_yBad->GetYaxis()->SetTitleOffset(1.4);
        TH2_radius_xRadius_yBad->GetXaxis()->SetNdivisions(505);

	TH1* TH1_F1 = new TH1F("TH1_F1","",50,0,.7);
	TH1_F1->SetTitle(";F_1 Score;");
	TH1_F1->SetLineWidth(2);

        // Create TTreeReader from this chain
        TTreeReader reader(SelectionTree);

        // Create TTreeReaderValues for the relevant branches
        TTreeReaderValue<double> radius_x(reader, "radius_x");
        TTreeReaderValue<double> radius_y(reader, "radius_y");
        TTreeReaderValue<double> z_1(reader, "z_1");
        TTreeReaderValue<double> z_2(reader, "z_2");
	TTreeReaderValue<double> mean_true_containment(reader, "mean_true_containment");
        TTreeReaderValue<double> min_shower_containment(reader, "min_shower_containment");
	TTreeReaderValue<double> min_shower_energy(reader, "min_shower_energy");
        TTreeReaderValue<double> min_shower_length(reader, "min_shower_length");
        TTreeReaderValue<double> max_shower_length(reader, "max_shower_length");
        TTreeReaderValue<double> max_shower_open_angle(reader, "max_shower_open_angle");

	TTreeReaderValue<double> NacceptedEvents(reader, "NacceptedEvents");
	TTreeReaderValue<double> truePositive(reader, "truePositive");
	TTreeReaderValue<double> falsePositive(reader, "falsePositive");
	TTreeReaderValue<double> trueNegative(reader, "trueNegative");
	TTreeReaderValue<double> falseNegative(reader, "falseNegative");

	double maximumPrecision = 0;
	double maximumF1 = 0;
	while(reader.Next())
	{
		double precision = *truePositive/(*truePositive+*falsePositive);
		double recall  = *truePositive/(*truePositive+*falseNegative);
		double F1 = 2 / (1/precision+1/recall);
		//if (precision < .96) 
		//{
		//	TH2_radius_xRadius_yBad->Fill(*radius_x, *radius_y);
		//	continue;
		//}
		if (precision > maximumPrecision) maximumPrecision = precision;
		if (F1 > maximumF1) maximumF1 = F1;
		TH2_radius_xRadius_y->Fill(*radius_x, *radius_y);
		TH2_precisionMinShowerContainment->Fill(*min_shower_containment, precision);
		TH2_recallMinShowerContainment->Fill(*min_shower_containment, recall);
		TH2_precisionRecall->Fill(precision, recall);
		TH1_F1->Fill(F1);	
	}

	cout<<"Maximum Precision: "<<maximumPrecision<<endl;
	cout<<"Maximum F1: "<<maximumF1<<endl;
	
	TCanvas *can = new TCanvas("can", "can", 800, 800);
	
	TH2_precisionRecall->Draw("COLZ");
	can->Print(Form("plots/%sTH2_Precision_vs_Recall.png", output_prefix.c_str()) );

	TH2_radius_xRadius_y->Draw("SCAT");
	TH2_radius_xRadius_yBad->Draw("SCAT same");
        //can->Print(Form("plots/%sTH2_Radius_x_vs_Radius_y.png", output_prefix.c_str()) );
	
	TH1_F1->Draw("hist");
	can->Print(Form("plots/%sTH1_F1score.png", output_prefix.c_str()) ); 

	can->cd();
	TPad *pad = new TPad("pad","",0,0,1,1);
	pad->Draw();
	pad->cd();
	TH2_precisionMinShowerContainment->Draw("SCAT");
	
	TPad *overlay = new TPad("overlay", "",0,0,1,1);
	overlay->SetFillStyle(4000);
	overlay->SetFrameFillStyle(4000);
	overlay->Draw();
	overlay->cd();	
	TH2_recallMinShowerContainment->Draw("SCAT Y+");
	//can->Print( Form("plots/%sTH2_MinShowerContainment_vs_precisionAndRecall.png", output_prefix.c_str()) ); 

	delete can;
	delete TH1_F1;
	delete TH2_precisionRecall;
	delete SelectionTree;
	delete f;
}


/*
 * Create 2D Histogram of reco and true shower containment for selected and rejected 
 * events using specified cuts
 */
void SelectionRegion(string output_prefix = "")
{
        // If output_prefix is not empty, append underscore after
        if (output_prefix.length() > 0)
                output_prefix.append("_");

        // Set global variables for ROOT
        gROOT->Reset();
        gStyle->SetOptStat(0);
        gStyle->SetLegendBorderSize(0);

        TFile *f = new TFile("CompressedPandoraLEEShower.root", "READ");
        TTree* MCsample = (TTree*)f->Get("ShowerAnalysisTree");

        // Create chain of ROOT file samples by reading in all files in nue.txt
        Long64_t N_events = MCsample->GetEntries();
        cout<<MCsample->GetName()<<" has "<<N_events<<" events"<<endl;

        // Create TTreeReader from this chain
        TTreeReader reader(MCsample);

        // Create TTreeReaderValues for the relevant branches
        // True shower energy (fiducial) obtained by looping through all daughters
        TTreeReaderValue<double> true_E(reader, "true_shower_E");
        TTreeReaderValue<double> true_E_fid(reader, "true_shower_E_fiducial");

        // Reconstructed shower start vertex
        TTreeReaderValue<double> start_x(reader, "shower_x");
        TTreeReaderValue<double> start_y(reader, "shower_y");
        TTreeReaderValue<double> start_z(reader, "shower_z");

        // Total energy in Y plane length and open angle
        // of shower closest to true vertex
        TTreeReaderValue<double> shower_energy(reader, "shower_energy");
        TTreeReaderValue<double> shower_length(reader, "shower_length");
        TTreeReaderValue<double> shower_open_angle(reader, "shower_open_angle");

        // True electron energy and shower containment
	TTreeReaderValue<double> true_primary_E(reader, "true_primary_e_energy");
        TTreeReaderValue<double> shower_containment(reader, "shower_containment");

	TH2* TH2_TrueRecoAgreementPass = new TH2F("TH2_TrueRecoAgreementPass","",50,0,1.001,50,0,1.001);
	TH2_TrueRecoAgreementPass->SetTitle("Selected Events;True Shower Containment;Reco Shower Containment");
	TH2_TrueRecoAgreementPass->GetYaxis()->SetTitleOffset(1.4);
	TH2_TrueRecoAgreementPass->SetMarkerStyle(20);
	TH2_TrueRecoAgreementPass->SetMarkerSize(.35);
	TH2_TrueRecoAgreementPass->SetMarkerColorAlpha(kBlue, 0.35);

	TH2* TH2_TrueRecoAgreementFail = new TH2F("TH2_TrueRecoAgreementFail","",50,0,1.001, 50,0,1.001);
        TH2_TrueRecoAgreementFail->SetTitle("Rejected Events;True Shower Containment;Reco Shower Containment");
	TH2_TrueRecoAgreementFail->GetYaxis()->SetTitleOffset(1.4);
	TH2_TrueRecoAgreementFail->SetMarkerStyle(20);
	TH2_TrueRecoAgreementFail->SetMarkerSize(.35);
        TH2_TrueRecoAgreementFail->SetMarkerColorAlpha(kRed, 0.35);

	TLine *bound1 = new TLine(0.05,0,1,0.95);
	TLine *bound2 = new TLine(0,0.05,0.95,1);
	bound1->SetLineColor(kBlack);
	bound2->SetLineColor(kBlack);
	bound1->SetLineWidth(2);
	bound2->SetLineWidth(2);
	
	
	// Cut parameters for Tight True Containment
	double trueContainmentThreshold = 0.95;
	double radius_x = 10;
	double radius_y = 25;
	double z_1 = 25;
        double z_2 = 870;
	double min_shower_containment = 0.9;
	double min_shower_energy = 0.25;
	double min_shower_length = 60;
	double max_shower_length = 225;
	double max_shower_open_angle = 0.5;
	
	/*
	// Cut parameters for Loose True Containment
        double trueContainmentThreshold = 0.8;
        double radius_x = 100;
        double radius_y = 85;
        double z_1 = 25;
        double z_2 = 950;
        double min_shower_containment = 0.7;
	*/

	double NacceptedEvents = 0;
	double truePositive = 0;
	double falsePositive = 0;
	double trueNegative = 0;
	double falseNegative = 0;

	while(reader.Next())
        {
        	// Skip if true shower differs from true electron energy by more than 10 MeV
                if (Abs(*true_primary_E-*true_E) > .01) continue;
                // Skip events when the shower begins in the dead region of the Y-plane
                if (*start_z >= 675. && *start_z <=775.) continue;
		
		NacceptedEvents++;

                double trueContainment = (*true_E_fid)/(*true_E);
                double Q = Power((*start_x-125.)/radius_x,2) + Power(*start_y/radius_y,2);
                if ( *start_z >= z_1 && 
	             *start_z <= z_2 && 
		     //Q <= 1. &&
		     //((*start_z >= 10 && *start_z <= 675) || (*start_z>=775 && *start_z <=986.8)) &&
		     *start_x >= 10 && *start_x <= 246.35 &&
		     Abs(*start_y) <= 96.5 && 
		     *shower_containment >= min_shower_containment )// &&
	   	     //*shower_energy >= min_shower_energy && 
		     //*shower_length >= min_shower_length && 
		     //*shower_length <= max_shower_length && 
		     //*shower_open_angle <= max_shower_open_angle )
                {
			//if (Abs(trueContainment-*shower_containment) <= trueRecoAgreement)
                	if (trueContainment >= trueContainmentThreshold) 
				truePositive++;
                        else falsePositive++;
			TH2_TrueRecoAgreementPass->Fill(trueContainment, *shower_containment);
		}
                else
                {
			//if (Abs(trueContainment-*shower_containment) <= trueRecoAgreement)
                        if (trueContainment < trueContainmentThreshold) 
				trueNegative++;
                        else falseNegative++;
                	TH2_TrueRecoAgreementFail->Fill(trueContainment, *shower_containment);
		}
        }

	/*
	cout<<"Total number of events: "<<NacceptedEvents<<endl;
	cout<<"Precision: "<<truePositive/(truePositive+falsePositive)*100.<<"%"<<endl;
	cout<<"Recall: "<<truePositive/(truePositive+falseNegative)*100.<<"%"<<endl;
	cout<<"Efficiency: "<<(truePositive+falsePositive)/NacceptedEvents*100.<<"%"<<endl;
	cout<<"Total number of selected events: "<<truePositive+falsePositive<<endl;
	cout<<"Average true containment of selected region: "<<TH1_TrueContainmentPass->GetMean()*100.<<"%"<<endl;
	cout<<"RMS of true containment of selected region: "<<TH1_TrueContainmentPass->GetRMS()*100.<<"%"<<endl;
	cout<<"Average true containment of rejected region: "<<TH1_TrueContainmentFail->GetMean()*100.<<"%"<<endl;
        cout<<"RMS of true containment of rejected region: "<<TH1_TrueContainmentFail->GetRMS()*100.<<"%"<<endl;
	*/

	TCanvas *can2 = new TCanvas("can","can",1700, 800);
	can2->Divide(2,1);
	can2->cd(1);
	gPad->SetRightMargin(0.15);
	TH2_TrueRecoAgreementPass->Draw("colz");
	bound1->Draw();
	bound2->Draw();
	can2->cd(2);
	gPad->SetRightMargin(0.15);
        TH2_TrueRecoAgreementFail->Draw("colz");
	bound1->Draw();
	bound2->Draw();
        can2->Print(Form("plots/%sTrueRecoAgreement.png", output_prefix.c_str()) );
}
