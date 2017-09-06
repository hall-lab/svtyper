develop:
	pip install -r requirements.txt --use-mirrors

test:
	python setup.py test

clean:
	find tests -name "*.pyc" -exec rm -v {} \;
	find scripts -name "*.pyc" -exec rm -v {} \;
	find svtyper -name "*.pyc" -exec rm -v {} \;
	if [ -e svtyper.pyc ]; then rm -v svtyper.pyc; fi;
	if [ -e .eggs ]; then rm -rfv .eggs; fi;
