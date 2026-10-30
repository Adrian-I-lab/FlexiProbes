# Security policy

## Supported versions

Only the latest release receives security fixes.

| Version | Supported |
|---|---|
| 1.0.x | Yes |
| < 1.0 | No |

## Reporting a vulnerability

Please report it privately, not as a public issue, using GitHub's private
reporting form:
https://github.com/Adrian-I-lab/FlexiProbes/security/advisories/new

Include the FlexiProbes version, your operating system, and the steps that
show the problem.

FlexiProbes is maintained by one person. I will do my best to reply promptly,
but please allow 2 to 4 weeks for a reply. If the report is accepted, the fix
ships in a new release and the advisory is published, crediting you unless
you ask not to be named. If it is declined, you will be told why.

## Scope

In scope: flaws in FlexiProbes's own code or packaging. Examples are an input
file (such as a FASTA or config file) that makes FlexiProbes run commands or
write outside the chosen output folder, or a generated HPC job script that
runs something the user did not ask for.

The bundled BLAST+, MAFFT and Qt are unmodified upstream builds. Report flaws
in them to their own projects, and tell me too so the bundled copy can be
updated.

Not a security issue: an incorrectly designed probe. Please use the
Issues template for reporting bugs.
