default: py check

.PHONY: py check

py:
	$(MAKE) --print-directory -C ..


# auto apply autopep8 and then
check:
	command -v autopep8 && autopep8 --in-place *.py --select=E261,W291,W293,W391 || echo "skipped autopep8"
	command -v flake8 && flake8 *.py
	command -v pycodestyle && for file in *.py ; do pycodestyle --ignore=E731,W503,W504,E129 $$file || exit; done || echo "skipped pycodestyle"
	@echo "done."

