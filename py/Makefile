all: build install

build: *.cpp *.h junkoda_cellularlib/*.py
	python3 setup.py config build_ext --inplace

install:
	pip install -e .

.PHONY: clean install check

clean:
	rm -rf build dist cellularlib.egg-info

check:
	$(MAKE) check --print-directory -C cellularlib
