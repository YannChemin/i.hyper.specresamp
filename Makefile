MODULE_TOPDIR = ../..

PGM = i.hyper.specresamp

include $(MODULE_TOPDIR)/include/Make/Script.make

default: script html

# Install the HTML manual
$(HTMLDIR)/$(PGM).html: $(PGM).html
	$(INSTALL_DATA) $(PGM).html $(HTMLDIR)/$(PGM).html

# Override htmldir target to ensure manual gets installed
htmldir: $(HTMLDIR)/$(PGM).html
