# Implementation specifications

Specifications are grouped by algorithm or source project, not by repository-wide migration phase. Future integrations belong in their own sibling directories rather than extending an existing algorithm's specifications.

- [Deen](deen/README.md)

Each group separates mathematical operators (`kernel-*.md`), host bindings (`plugin.md`), acceptance criteria (`acceptance.md`), and reference evidence (`reference-mapping.md`). Extract shared specifications only when multiple algorithms actually share the same semantics; no common layer is reserved in advance.

These documents define behavioral requirements independently of implementation status. Specification IDs do not depend on the neo-smo namespace, development phases, or source file locations. Other projects may reuse the algorithm contracts with their own host bindings.
