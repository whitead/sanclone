# SanClone


## Installation

```sh
pip install -e .
```

## Testing

```sh
pip install -r dev-requirements.txt
pytest
```

## Developing

To contribute, make sure to run the pre-commit hooks before committing.


### First time setup
```sh
# install pre-commit
pip install -r dev-requirements.txt
pre-commit install
pre-commit run --all-files
```

After that, it will run automatically
