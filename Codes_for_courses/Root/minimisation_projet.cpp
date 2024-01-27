#include <iostream>

//On crée la fonction construct() qui construit
//Les histogrammes associés au bruit et au signal
//et à la somme des deux:
void construct(){

    //Variables:
    double amplitude_bruit = 1;                                 //Amplitude du fond  
    double slope_bruit = 1.5;                                   //x dans l'exponentielle du bruit
    double n_signal = 100;                                      //n_total*k;//Nombre de répétitions du signal
    double n_bruit = 1000;                                      //n_total-n_signal;
    double sigma = 0.175;                                       //Déviation std de la gaussienne du signal
    double mu = 1.5;                                            //Moyenne de la gaussienne du signal
    double nbins = 100;                                         //Nombre de bins
    double amplitude_signal = 1/sqrt(2*TMath::Pi()*sigma);      //Amplitude du signal
    
    TF1* bruit = new TF1("bruit","[0]*exp(-[1]*x)",0,10);       //On crée une fonction de bruit exp décroissante
    bruit->SetParameters(amplitude_bruit,slope_bruit);          //On introduit les paramètres

    TF1* signal = new TF1("signal","gaus",0,10);                //On crée le signal gaussien
    signal->SetParameters(amplitude_signal,mu,sigma);           //On ajuste les paramètres

    TH1F* htotal = new TH1F("total","total",nbins,0,10);        //On initialise l'histogramme total du signal
    htotal->FillRandom("bruit",n_bruit);                        //On rempli avec les répétitions
    htotal->FillRandom("signal",n_signal);                      //On rempli avec les répétitions
    //htotal->Draw("e");
    
    //Fichier de sortie de l'Histograme total:
    TFile fout("total.root","recreate");
	htotal->Write();
	fout.Close();
}

//On crée la fonction likelihood_imp() qui implémente
//la méthode de minimum de vraissemblance sur les
//signaux

double likelihood_imp(double k){ 
/************************ Initialisation ****************************/

    //On lit les fichiers des histogrammes:
    TFile ftotal("total.root");//Histogramme total
    
    //On garde les histogrammes dans des variables
    TH1F* htotal = (TH1F*)ftotal.Get("total");                  //Histogramme total

    //On déclare les variables nécessaires pour l'analyse
    
    int nbins = htotal->GetNbinsX();                            //Le nombre de bins
    
    //Quantités de répétitions dans chaque bin:
    double bintotal[nbins];                                     //Dans l'histogramme total
    //Centre des bins:
    double center[nbins];                                       //Pour l'histogramme total
    double dist[nbins];                                         //Distribution de poisson associée à k
    double l=0;                                                 //-log du likelihood pour chaque k
    double mean[nbins];                                         //Espérance de la distribution de Poisson
    double sigma = 0.175;                                       //Déviation std de la gaussienne du signal
    double amplitude_signal = 1/sqrt(2*TMath::Pi()*sigma);      //Amplitude du signal
    double mu = 1.5;                                            //Moyenne de la gaussienne du signal
    double Ntot = htotal->GetEntries();                         //Nombre total d'événements
    double largeur = htotal->GetBinWidth(1);                    //Largeur des bins
    double slope_bruit = 1.5;                                   //Comme son nom l'indique
    TF1* model = new TF1("model","[0]*([1]*gaus(2)+(1-[1])*[5]*exp(-[6]*x))",0,10);//On crée le signal gaussien
    model->SetParameters(Ntot*largeur,k,amplitude_signal,mu,sigma,1,slope_bruit);//On ajuste les paramètres

/************************* Application de la méthode ************************************/
    for(int i=0;i<nbins;i++){                                   //on réalise une boucle sur les bins
        bintotal[i] = htotal->GetBinContent(i);                 //Nombre d'évènements totaux par bin
        center[i] = htotal->GetBinCenter(i);                    //Centre des bins
        mean[i] = model->Eval(center[i]);                       //Espèrance des évènements
        dist[i] = TMath::Poisson(bintotal[i],mean[i]);          //On calcule la distribution de Poisson
        l += -2*log(dist[i]);                                   //On calcule la fonction recherchée
    }
    return l;
}

int minimisation(){
    double kmax = 0.12;                                         //Borne supérieure de k
    double kmin = 0.07;                                         //Borne inférieure de k
    int nrep = 100;                                             //nombre de bins
    double deltak = (kmax-kmin)/nrep;                           //pas sur k
    double k[nrep];                                             //vecteur k
    double likelihood[nrep];                                    //-2log(likelihood)
    int Ndata = 1000;                                           //Nombre de lots de données générées
    double pull[Ndata];                                         //Vecteur de pulls
    TH1F* hpull = new TH1F("pull","pull",50,5,15);              //Histogramme des pulls
    for(int n=0;n<Ndata;n++){                                   //On boucle sur le nombre de lots générés
        construct();                                            //On génére les données
        for(int i=0;i<nrep;i++){                                //On boucle sur les bins pour calculer les likelihood
            k[i]=kmin+deltak*i;      
            likelihood[i]=likelihood_imp(k[i]);    
        }

        
        TGraph *gr = new TGraph(nrep,k,likelihood);             //On crée le graph
        //TCanvas *c1 = new TCanvas("c1","Graph Draw Options",200,10,600,400);
        gr->Fit("pol2");                                        //On fit
        //gr->Draw("A*");
        //gr->Draw("C");
        //On obtient les parmètres du fit:
        TF1* fitted = gr->GetFunction("pol2");
        double a = fitted->GetParameter(2);
        double b = fitted->GetParameter(1);
        double c = fitted->GetParameter(0);
        double minim = -b/(2*a);                                //Le minimum du likelihood
        double delta = 1/sqrt(2*a);                             //La déviation standard de l'estimateur
        double deltamax = minim+delta;                          //Borne supérieure de l'intervalle de confiance
        double deltamin = minim-delta;                          //Borne inférieure de l'intervalle de confiance

        pull[n] = (minim-100/1100)/(delta);                     //On calcule les pulls
        hpull->Fill(pull[n]);                                   //On crée l'histogramme des pull
    }
    //On trace l'histogramme:
    TCanvas *c1 = new TCanvas("c1","Graph Draw Options",200,10,600,400);
    hpull->Fit("gaus");
    hpull->Draw();
    return 0;

}