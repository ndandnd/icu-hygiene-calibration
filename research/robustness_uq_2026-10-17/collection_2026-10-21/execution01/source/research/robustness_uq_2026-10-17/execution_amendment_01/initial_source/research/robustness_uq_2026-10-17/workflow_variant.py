"""Isolated vacancy-cleaning/empty-isolation-bed robustness model.

The maintained kernel and earlier frozen studies are not modified. This module
loads the hash-guarded passive response instrumentation and subclasses only the
discharge and relocation workflow. It is an alternative illustrative workflow,
not a claim that these rules reproduce a particular ICU.

Use ``ReferenceModel`` for the unchanged response kernel, or pass
``WorkflowRobustnessModel`` as ``model_factory`` to the response pair runner.
Additional action counters are carried by its existing snapshot/difference
mechanism. No additional random draws are made by these workflow operations.
"""

from collections import Counter
import hashlib
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
ABM = HERE.parents[1]
INSTRUMENTATION_PATH = ABM / "research/full_icu_response_2026-10-08/instrumented_model.py"
INSTRUMENTATION_SHA256 = "5dcbe23b392694c1baffa7cd2475588ab12d8a3333778c17e73e0b5e6131888b"

if hashlib.sha256(INSTRUMENTATION_PATH.read_bytes()).hexdigest() != INSTRUMENTATION_SHA256:
    raise RuntimeError("Frozen response instrumentation changed")

_MODULE_NAME = "_robustness_uq_response_instrumentation"
if _MODULE_NAME in sys.modules:
    instrumentation = sys.modules[_MODULE_NAME]
    if Path(instrumentation.__file__).resolve() != INSTRUMENTATION_PATH.resolve():
        raise RuntimeError("Response instrumentation module identity collision")
else:
    _spec = importlib.util.spec_from_file_location(_MODULE_NAME, INSTRUMENTATION_PATH)
    instrumentation = importlib.util.module_from_spec(_spec)
    sys.modules[_MODULE_NAME] = instrumentation
    _spec.loader.exec_module(instrumentation)

ReferenceModel = instrumentation.InstrumentedModel


class WorkflowRobustnessModel(ReferenceModel):
    """Clear vacated bedside sites and relocate only to empty isolation beds.

    Cleaning is a deterministic full site clear immediately after discharge or
    relocation. It clears the site's lineage but leaves its scheduled-clean
    countdown unchanged. The destination is not additionally cleaned at move-in:
    any contamination acquired while it was vacant is retained. Fixed bed order,
    screening, admission, staffing, hand states, and all other dynamics remain
    inherited. Failed relocation attempts are counted once per eligible patient
    per invocation (normally once per tick), not as distinct patients.
    """

    workflow_name = "clean_vacancy_empty_isolation"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        sites_by_zone = {
            self.care_zone_by_position[site.pos]: site
            for site in self.response_sites
        }
        self.response_site_by_bed_position = {
            bed.pos: sites_by_zone[bed.unique_id]
            for bed in self.shared_beds + self.isol_beds
        }
        if len(sites_by_zone) != 30 or len(self.response_site_by_bed_position) != 30:
            raise RuntimeError("Workflow requires one fixed environmental site per bed")
        self.response_counts["workflow_events"] = Counter()

    def _clear_vacated_site(self, bed_position, reason):
        if reason not in ("discharge", "relocation"):
            raise ValueError("Unknown vacancy-cleaning reason")
        if any(a.isPatient for a in self.grid.get_cell_list_contents([bed_position])):
            raise RuntimeError("Cannot clean a site as vacant while its bed is occupied")
        site = self.response_site_by_bed_position[bed_position]
        was_dirty = bool(site.colonized)
        site.colonized = False
        site._response_environment_origin = None
        self.response_environment_count -= int(was_dirty)
        events = self.response_counts["environment_events"]
        events[reason + "_clear"] += 1
        if was_dirty:
            events["dirty_" + reason + "_clear"] += 1

    def _remove_discharged(self):
        vacated_positions = [patient.pos for patient in self.discharged]
        super()._remove_discharged()
        for position in vacated_positions:
            self._clear_vacated_site(position, "discharge")

    def _move_patient(self, patient, destination):
        origin = patient.pos
        if destination == origin:
            return
        if destination not in self.response_site_by_bed_position:
            raise ValueError("Patient destination must be a fixed bed")
        if any(a.isPatient for a in self.grid.get_cell_list_contents([destination])):
            self.response_counts["workflow_events"]["occupied_destination_rejected"] += 1
            raise ValueError("Workflow relocation requires an empty destination bed")
        super()._move_patient(patient, destination)
        self._clear_vacated_site(origin, "relocation")
        self.response_counts["workflow_events"]["isolation_relocation"] += 1

    def _relocate_detected_patients(self):
        for bed in self.shared_beds:
            occupants = [
                a for a in self.grid.get_cell_list_contents([bed.pos])
                if a.isPatient and a.positive
            ]
            if not occupants:
                continue
            patient = occupants[0]
            for isolated_bed in self.isol_beds:
                isolated_bed.checkFilled()
                if isolated_bed.filled:
                    continue
                self._move_patient(patient, isolated_bed.pos)
                bed.checkFilled()
                isolated_bed.checkFilled()
                break
            else:
                self.response_counts["workflow_events"]["isolation_deferred_no_empty_bed"] += 1

    def assert_response_invariants(self):
        super().assert_response_invariants()
        environment = self.response_counts["environment_events"]
        workflow = self.response_counts["workflow_events"]
        assert environment["discharge_clear"] == self.discharged_patients
        assert environment["relocation_clear"] == workflow["isolation_relocation"]
        for reason in ("discharge", "relocation"):
            assert environment["dirty_" + reason + "_clear"] <= environment[reason + "_clear"]


__all__ = ["ReferenceModel", "WorkflowRobustnessModel"]
