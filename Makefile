.PHONY: default clean check push

default:
	pip install -e py

clean:
	$(MAKE) clean --print-directory -C py

push:
	sh push.sh
