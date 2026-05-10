ifeq ($(DEBIAN_BUILD),1)
PREFIX   ?= /usr/local
BINDIR    = $(DESTDIR)$(PREFIX)/bin
PGM_PY   = $(PGM).py

$(PGM): $(PGM_PY)
	cp $< $@ && chmod 755 $@

install: $(PGM)
	install -d $(BINDIR)
	install -m 755 $(PGM) $(BINDIR)/$(PGM)

clean:
	rm -f $(PGM)

.PHONY: install clean

else
MODULE_TOPDIR = ../..

PGM = i.hyper.specresamp

include $(MODULE_TOPDIR)/include/Make/Script.make

default: script html

# Install the HTML manual
$(HTMLDIR)/$(PGM).html: $(PGM).html
	$(INSTALL_DATA) $(PGM).html $(HTMLDIR)/$(PGM).html

# Override htmldir target to ensure manual gets installed
htmldir: $(HTMLDIR)/$(PGM).html
endif
