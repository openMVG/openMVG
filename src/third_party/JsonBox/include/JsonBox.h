#ifndef JB_JSON_BOX_H
#define JB_JSON_BOX_H

/**
 * @mainpage JsonBox
 * @section format_description What is JSON?
 * [JSON](http://json.org/) (JavaScript Object Notation) is a lightweight data-interchange format.
 * @section project_description What is Json Box?
 * Json Box is a C++ library used to read and write JSON with ease and speed.
 *
 * Things it does:
 * * Follows the standards established on [http://json.org/](http://json.org/)
 * * Read and write JSON in UTF-8
 * * Uses the STL streams for input and output
 * * Generated JSON can be indented and pretty or compact and hard-to-read
 * * Does not crash when the JSON input contains errors, it simply tries to interpret as much as it can
 *
 * Things it does not do:
 * * Read JSON in UTF-16 or UTF-32
 * * Keep the order of the members in objects (the standard doesn't require keeping the order)
 * * Write useful error messages when the JSON input contains errors
 *
 * The library wasn't designed with multi-threading in mind.
 * @see JsonBox
 */

#include <JsonBox/Value.h>
#include <JsonBox/Array.h>
#include <JsonBox/Object.h>

#endif
