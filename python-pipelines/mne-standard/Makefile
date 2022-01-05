
help:
	@echo "Please use \`make <target>' where <target> is one of"
	@echo "  ica           use ICA to remove eye blinks and heart beat"
	@echo "  epochs        create epochs and evoked responses"
	@echo "  all           all processing "


SUBJ := $(shell cat subjects.txt)

all: ica epochs forward inverse

tsss-ses1:
	for i in $(SUBJ); do \
		python tsss.py $$i -ses 01 -task envspeech -run 01; \
		python tsss.py $$i -ses 01 -task restEO -run 01; \
		python tsss.py $$i -ses 01 -task restEC -run 01; \
	done

tsss-ses2:
	for i in $(SUBJ); do \
		python tsss.py $$i -ses 02 -task envspeech -run 01; \
		python tsss.py $$i -ses 02 -task restEO -run 01; \
		python tsss.py $$i -ses 02 -task restEO -run 02; \
	done

ica:
	for i in $(SUBJ); do \
		python ica.py $$i;\
	done

epochs:
	for i in $(SUBJ); do \
        python epochs.py $$i -ses 01 -task envspeech;\
        python epochs.py $$i -ses 02 -task envspeech;\
	done

forward: 
	for i in $(SUBJ); do \
		python forward.py $$i;\
	done

noise_cov: 
	for i in $(SUBJ); do \
		python noise_cov.py $$i -ses 01;\
		python noise_cov.py $$i -ses 02;\
	done	

inverse:
	for i in $(SUBJ); do \
		python inverse.py $$i -ses 01 -task envspeech;\
		python inverse.py $$i -ses 02 -task envspeech;\
	done

#.PHONY clean
#clean:
#	python clean.py

