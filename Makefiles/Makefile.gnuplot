################################################################################
# Compile rules - Gnuplot
################################################################################
# Gnuplot tools processing rules
#

$(GNUPLOT_ROOT)/Makefile:
	cd $(GNUPLOT_ROOT); ./configure --without-readline --without-lisp-files --without-tutorial --without-row-help --without-x --disable-wxwidgets --with-png=$(realpath $(LIBPNG_ROOT)) --with-zlib=$(realpath $(LIBZ_ROOT)) --with-gd=$(realpath $(LIBGD_ROOT)) --with-jpeg --with-freetype $(MINGW_SPECIFIC_FLAGS); cd $(ROOT)
	#cd $(GNUPLOT_ROOT); ./configure --with-jpeg --with-png --with-freetype ; cd $(ROOT)

$(GNUPLOT_ROOT)/src/gnuplot: $(GNUPLOT_ROOT)/Makefile
	$(MAKE) -C $(GNUPLOT_ROOT)/src LDFLAGS="$(LDFLAGS2)" CFLAGS="$(CPPFLAGS) -fPIC"
	# $(MAKE) -C $(GNUPLOT_ROOT)/src LDFLAGS="$(LDFLAGS) -lpthread -static-libgcc -lexpat" CFLAGS="$(CPPFLAGS) -fPIC"

$(PREFIX)/$(SUFFIX)/gnuplot: $(GNUPLOT_ROOT)/src/gnuplot
	cp $(GNUPLOT_ROOT)/src/gnuplot $(PREFIX)/$(SUFFIX)

$(BIN_DIR)/gnuplot: # $(PREFIX)/$(SUFFIX)/gnuplot
	cp $(PREFIX)/$(SUFFIX)/gnuplot $(BIN_DIR)

gnuplotExe: $(BIN_DIR)/gnuplot

################################################################################
