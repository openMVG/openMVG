<!-- -*- markdown -*- -->
# Json Box
[JSON](http://json.org/) (JavaScript Object Notation) is a lightweight data-interchange format.

Json Box is a C++ library used to read and write JSON with ease and speed.

Things it does:
* Follows the standards established on [http://json.org/](http://json.org/)
* Read and write JSON in UTF-8
* Uses the STL streams for input and output
* Generated JSON can be indented and pretty or compact and hard-to-read
* Does not crash when the JSON input contains errors, it simply tries to interpret as much as it can

Things it does not do:
* Read JSON in UTF-16 or UTF-32
* Keep the order of the members in objects (the standard doesn't require keeping the order)
* Write useful error messages when the JSON input contains errors

The library wasn't designed with multi-threading in mind.

The class reference can be found [here](http://anhero.github.com/JsonBox).

