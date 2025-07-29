sphinx-quickstart docs
sphinx-apidoc -o docs/source/GenoKit GenoKit
sphinx-build -b html docs/source/ docs/build/html
cd docs
sphinx-build -b html source/ build/
make html
