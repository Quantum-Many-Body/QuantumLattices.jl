# Code style guide

Basically, the code style is compatible with `Julia.Base`, which follows much of the rules that is also popular in other languages.

## Whitespaces

* Always surround these binary operators with a single space on either side: assignment (=), augmented assignment (+=, -= etc.), comparisons (==, <, >, !=, <>, <=, >=, in, not in, is, is not), Booleans (and, or, not)

* No white spaces immediately inside parentheses, brackets or braces.
  ```
  Yes: sum(1)
  No:  sum( 1), sum(1 ), sum( 1 )
  ```

## Naming Conventions