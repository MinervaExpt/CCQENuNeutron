CXX := g++
ROOTFLAGS := `root-config --cflags --glibs --libs`
CXXFLAGS := -g -Wall -fPIC 
# Make sure your ROOTFLAGS includes -lTreePlayer, which is where the proxy magic lives
INCLUDE += -I ./ -I $(UNFOLDUTILSROOT) -I $(PLOTUTILSROOT)/ -I $(CCQENU_ANA)/ -I$(CCQENU_INCLUDE_PATH) -I$(MCREWEIGHTROOT) -I $(MY_CCQENU) 
LDLIBS += -L./ -L$(UNFOLDUTILSROOT)/$(CMTCONFIG) -lUnfoldUtils -L$(PLOTUTILSROOT)/$(CMTCONFIG) -lplotutils -lCintex -lMinuit -lGenVector -L$(MCREWEIGHTROOT)/$(CMTCONFIG) -lMCReweight -L$(CCQENU_ANA)/lib -lCCQENU
# Add any more executables you want here. They should have a single
# .cxx file with the same name as the executable you want

#BINARIES = plot_antinu_xsecs plot_nu_xsecs plot_nu_xsecs2 error_plot_grid plot_nu_xsecs_models plot_antinu_xsecs2 plot_antinu_xsecs_models plot_nu_xsecs_models_datacrosssection plot_nu_varied_iterations plot_nu_varied_iterations_ratio plot_nu_vtx_energy plot_nu_vtx_energy_evt_type error_plot_grid_vtx_energy plot_nu_eventrate plot_nu_vtx_energy_qescale error_plot_grid_3D plot_nu_xsec3 plot_nu_xsec3_models plot_nu_xsec3_compact plot_nu_bkg plot_nu_transverse plot_nu_recoil_nuescatter error_plot_grid_transverse  plot_nu_eventrate_3D plot_nu_eventrate_inclusive plot_nu_xsec2_inclusive plot_nu_eff_3D plot_nu_xsecs2_ME plot_nu_xsecs2_ME_NS error_plot_grid_ME error_plot_grid_ME_NS_truncated error_plot_grid_ME_NS plot_nu_xsecs_models_MnvGENIE_ME  plot_nu_eventrate_ME plot_nu_eventrate_pion
#BINARIES = plot_nu_eventrate plot_nu_eventrate_neutron plot_nu_eventrate_neutron_proton_reweights plot_nu_eventrate_neutron_only_reweights plot_nu_eventrate_2D_neutron_proton_reweights plot_nu_eventrate_regions_neutron_proton_reweights 
BINARIES = plot_nu_eventrate1 plot_nu_eventrate2 plot_nu_eventrate3 plot_nu_eventrate4 plot_nu_eventrate5 plot_nu_eventrate_regions_neutron_proton_2D plot_nu_eventrate_regions_neutron_proton_eff_mig_xs plot_pn_eventrate_regions_neutron_proton_func plot_nu_diagnostic_plots plot_nu_eventrate_test  FATest
#plot_nu_xsecs_models_ME 

all : $(BINARIES) libGridCanvas.so

libGridCanvas.so: GridCanvas.o GridCanvas_Dict.o
	$(CXX) $(ROOTFLAGS) -shared -o $@ $^

%.o: %.cxx plot.h
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS)  -c $< -o $@ 


%.o: %.cxx plot_ME.h
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS)  -c $< -o $@ 

GridCanvas_Dict.cxx: GridCanvas.h LinkDef.h
	rootcint -f $@ -c $(CXXFLAGS) -p $^

# All of the binaries have the same format, so use a "static pattern
# rule". Each binary "foo" depends on "foo.o" and we build it with the
# recipe given ($@ will be the name of the binary)
$(BINARIES) : %: %.o libGridCanvas.so
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -lGridCanvas -o $@ $^

clean:
	rm -f $(BINARIES) *.o *.so GridCanvas_Dict.*

.PHONY: clean
