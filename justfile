# Gate run by the land queue before merging.
test-gate:
    make test
    make compile
    make checkdoc
