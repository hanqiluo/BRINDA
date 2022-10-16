# BRINDA 0.1.0

# BRINDA 0.1.1
- Resolve the problems of global variables;
- Change "=" to "<-" in some lines

# BRINDA 0.1.2
- Add more description
- Use "message()" to print notifications instead of "print()".
- Organize the notifications to make it flow better.
- Change all(is.na(agp)) or all(is.na(crp)) to exists("agp", dataset), or 
exists("crp", dataset)

# BRINDA 0.1.3
- Reduce the title to 65 characters
- Add more description and DOI

# BRINDA 0.1.4
- To adjust for inflammation in Zinc, correlation should < -0.1 instead of -0.2
- For zinc: applying adjustment for AGP when only AGP is available (instead of adjustment for AGP and CRP)
- Add more examples to the readme

# BRINDA 0.1.5
- Change the condition of adjusting for sTfR. sTfR is adjusted by AGP only in WRA and PSC population groups. sTfR can be adjusted by AGP or CRP in other/manual population groups.
- 10/16/2022：moved URL: https://brinda-nutrition.org/publications/ to
https://www.brinda-nutrition.org/publications)

