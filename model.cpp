#include <morph/HexGrid.h>
#include <morph/RD_Base.h>
#include <morph/Config.h>

#include <morph/HdfData.h>
#include <morph/Config.h>
#include <morph/ColourMap.h>
#include <morph/Visual.h>
#include <morph/HexGridVisual.h>
#include <morph/GraphVisual.h>
#include <morph/VisualDataModel.h>
#include <morph/Scale.h>
#include <morph/Vector.h>

typedef morph::VisualDataModel<FLT>* VdmPtr;

using morph::Config;
using morph::Tools;
using morph::ColourMap;


// DERIVED CLASS THAT IMPLEMENTS A SELF-PROJECTION - ADD FUNCTIONS TO IT IF YOU WANT!
template <class Flt>
class RD_Sheet : public morph::RD_Base<Flt>
{
public:

    float dt;
    std::vector<unsigned int> counts;                   // number of connections in connection field for each unit
    std::vector<std::vector<unsigned int> > neighbourID;      // identity of conneted units on the source sheet
    std::vector<std::vector<Flt> > weights;             // connection weights
    alignas(alignof(std::vector<Flt>)) std::vector<Flt> X;
    alignas(alignof(std::vector<Flt>)) std::vector<Flt> Xprev;
    std::vector<Flt> Xpos, Ypos;

    virtual void init (void) {
        this->stepCount = 0;
        this->zero_vector_variable (this->X);
    }

    virtual void allocate (void) {
        morph::RD_Base<Flt>::allocate();
        this->resize_vector_variable (this->X);
        this->resize_vector_variable (this->Xprev);

        for(int i=0; i<this->hg->vhexen.size();i++){
            Xpos.push_back(this->hg->vhexen[i]->x);
            Ypos.push_back(this->hg->vhexen[i]->y);
        }

    }

    void setSelfProjection(void){

        int N = this->hg->vhexen.size();
        counts.resize(N,0);
        neighbourID.resize(N);
        weights.resize(N);

        // initialize connections for each unit
        #pragma omp parallel for
        for(unsigned int i=0;i<N;i++){
            for(unsigned int j=0;j<N;j++){
                //Flt dx = (this->hg->vhexen[j]->x - this->hg->vhexen[i]->x);
                //Flt dy = (this->hg->vhexen[j]->y - this->hg->vhexen[i]->y);
                Flt dr = (this->hg->vhexen[j]->ri - this->hg->vhexen[i]->ri);
                Flt dg = (this->hg->vhexen[j]->gi - this->hg->vhexen[i]->gi);
                Flt db = (this->hg->vhexen[j]->bi - this->hg->vhexen[i]->bi);

                if(i==j){                           //connection from hex unit to itself
                    counts[i]++;
                    neighbourID[i].push_back(j);
                    weights[i].push_back(42.0);
                }
                if((abs(dr)+abs(dg)+abs(db)==1)){   //distance of 1 in hexgrid coords
                    counts[i]++;
                    neighbourID[i].push_back(j);
                    weights[i].push_back(-10.0);
                }
                if((abs(dr)+abs(dg)+abs(db)==2)){   //distance of 2 in hexgrid coords
                    counts[i]++;
                    neighbourID[i].push_back(j);
                    weights[i].push_back(1.5);      //note some of these should be 1 and some should be 2
                }

            }
        }

    }

    void step(void){
        Xprev = X;
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            X[hi] = 0.;
            for(unsigned int j=0;j<counts[hi];j++){
                X[hi] += Xprev[neighbourID[hi][j]]*weights[hi][j];
            }
            X[hi] *= dt;
        }
    }

    void noise_X (float gain) {
        morph::RandUniform<Flt> rng;
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->X[hi] = rng.get() * gain;
        }
    }

    ~RD_Sheet(){

    }


};




int main(int argc, char **argv){

    // GENERAL SETUP
    if (argc < 3) { std::cerr << "\nUsage: ./build/model configfile logdir seed\n\n"; return -1; }
    std::srand(std::stoi(argv[3]));       // set seed
    std::string paramsfile (argv[1]);
    Config conf(paramsfile);
    if (!conf.ready) { std::cerr << "Error setting up JSON config: " << conf.emsg << std::endl; return 1; }
    std::string logpath = argv[2];
    std::ofstream logfile;
    morph::Tools::createDir (logpath);
    { std::stringstream ss; ss << logpath << "/log.txt"; logfile.open(ss.str());}
    logfile<<"Hello."<<std::endl;

    // READ SETUP INFO FROM CONFIG FILE
    const unsigned int plotevery = conf.getUInt ("plotevery", 1);
    const bool saveplots = conf.getBool ("saveplots", false);
    unsigned int framecount = 0;
    const unsigned int win_height = conf.getUInt ("win_height", 400);
    const unsigned int win_width = conf.getUInt ("win_width", win_height);
    const unsigned int steps = conf.getUInt ("steps", 1000);

    // READ SIMULATION INFO FROM CONFIG FILE
    float noiseGain = conf.getFloat ("noiseGain", 0.1);
    float dt = conf.getFloat ("dt", 0.001);

    // INITIALIZE MODEL SHEET
    RD_Sheet<FLT> J;
    J.svgpath = conf.getString ("svgpath", "boundaries/trialmod.svg");
    J.init();
    J.allocate();
    J.setSelfProjection();
    J.noise_X(noiseGain);
    J.dt = dt;

    // SETUP PLOTTING
    std::chrono::steady_clock::time_point lastrender = std::chrono::steady_clock::now();
    morph::Visual v1 (win_width, win_height, "model");
    v1.backgroundWhite();
    v1.sceneLocked = conf.getBool ("sceneLocked", false);
    v1.scenetrans_stepsize = 0.1;
    v1.fov = 15;

    std::vector<unsigned int> grids(1);
    float txtoff = -0.6f;

    // ADD PLOTS TO SCENE

    morph::Scale<FLT> zscale; zscale.setParams (0.0f, 0.0f);
    morph::Scale<FLT> cscale; cscale.do_autoscale = true;

    morph::HexGridVisual<FLT> hgv (v1.shaderprog,v1.tshaderprog, J.hg,std::array<float,3>{0.0f,0.0f,0.0f}, &(J.X),zscale,cscale,morph::ColourMapType::Jet);
    grids[0] = v1.addVisualModel (&hgv);
    v1.getVisualModel (grids[0])->addLabel ("Hello Jay!!", {-0.5f, txtoff, 0.0f},
    morph::colour::black, morph::VisualFont::VeraSerif, 0.05, 56);
    hgv.hexVisMode = morph::HexVisMode::Triangles;


    // STEP THROUGH THE MODEL
    for(int t=0;t<steps; t++){

        // integrate dynamics
        J.step();

        //update plot
        if(t%plotevery==0){
            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids[0]);
            avm->updateData (&(J.X));
            avm->clearAutoscaleColour();
        }

        // get user input
        std::chrono::steady_clock::duration sincerender = std::chrono::steady_clock::now() - lastrender;
        if (std::chrono::duration_cast<std::chrono::milliseconds>(sincerender).count() > 17) {
            glfwPollEvents();
            v1.render();
            lastrender = std::chrono::steady_clock::now();
        }

        // save plot
        if(saveplots){
            std::stringstream ff;
            ff << logpath << "/model_";
            ff << std::setw(5) << std::setfill('0') << framecount;
            ff << ".png";
            v1.saveImage (ff.str());
            framecount++;
        }

    }


    // save data out to file
    {
        std::stringstream fname;
        fname << logpath << "/data.h5";
        morph::HdfData data(fname.str());
        std::stringstream ss;
        ss.str("");
        ss.clear();
        ss<<"/x";
        data.add_contained_vals (ss.str().c_str(), J.Xpos);
        ss.str("");
        ss.clear();
        ss<<"/y";
        data.add_contained_vals (ss.str().c_str(), J.Ypos);
        ss.str("");
        ss.clear();
        ss<<"/X";
        data.add_contained_vals (ss.str().c_str(), J.X);
    }


    std::cout<<"Finished."<<std::endl;
    return 0;
}
