SUBDIRS := background_tuning make_hists $(CCQENU_ANA)/src

all: $(SUBDIRS)

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done

background_tuning: $(CCQENU_ANA)/src
make_hists: $(CCQENU_ANA)/src
make_plots: $(CCQENU_ANA)/src
utils: $(CCQENU_ANA)/src

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: all clean $(SUBDIRS)

