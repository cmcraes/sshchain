# SSH Model Code

Experimental code for calculation of Lohschmitt echo 
in Su-Schrieffer-Heeger chains 

Since it keeps raining, I write some more documentation.

## Installation

To install this package,

```{bash}
cd /where/ever/you/put/it/sshcain
cmake .
make
```
To build the documentation use `doxygen` or `make doc`.

## Dependencies

### boost <= 1.53
is needed mainly for the `boost::multiprecision` wrappers

### GMP (which version??)
is the backend used by boost. Probably, the version shipped with mpack should work.

### mpalapck/mpack
You can build `mplapck` from scratch, but it is a pin in the neck. 
I have a binary package which is way to big for github, but I can gibe it upon request.

### GOMP
for reasons

