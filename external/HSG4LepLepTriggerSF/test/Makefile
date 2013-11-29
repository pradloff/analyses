all:
	g++ -o test test.cxx `root-config --cflags --libs` -L../StandAlone -lHSG4LepLepTriggerSF -I../HSG4LepLepTriggerSF
	g++ -o new_test new_test.cxx `root-config --cflags --libs` -L../StandAlone -lHSG4LepLepTriggerSF -I../HSG4LepLepTriggerSF

clean:
	-rm -f test new_test


