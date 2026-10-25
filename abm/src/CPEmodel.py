"""Corrected, site-specific CPE carriage model.

Historical experiments and calibration are archived, not validated by this code.
See README.md for intentional behavior changes and unresolved assumptions.
"""
import math
import numbers

import numpy as np
from mesa import Model
from mesa.datacollection import DataCollector
from mesa.space import MultiGrid
from mesa.time import BaseScheduler

from agents import Bed, Dr, Goo, IsolatedBed, Nurse, Patient, XrayDr, ICUA, ICUB

width = 32
height = 11


def getTotPatients(model):
    return len(model.current_patients)


def getNumSick(model):
    return sum(a.colonized for a in model.current_patients)


def getNumColonized(model):
    return sum(a.isHCW and a.colonized for a in model.schedule.agents)


def getCumulNumSick(model):
    return model.cumul_sick_patients


def getHCWInfec(model):
    """Legacy reporter name: this counts acquired colonization, not infection."""
    return model.acquired_colonizations


def getCumul(model):
    return model.cumul_patients


paths = np.concatenate((np.array_split(ICUA, 5), np.array_split(ICUB, 5)))


class CPE_Model(Model):
    """A 30-bed, ten-nurse illustrative carriage model with owned RNGs.

    Synthetic LOS defaults (5 and 14 days) are scenario assumptions, not clinical
    estimates. Historical LOS resampling is available only by explicit opt-in.
    Each tick represents 1/2592 day. Census exposure is counted at tick start.

    Hygiene completion has marginal probability given by each worker's wash
    rate, lag-one correlation ``hygiene_autocorrelation``, and conditional
    binary clearing probability ``hygiene_efficacy``. ``separate`` provides a
    distinct before and after opportunity per contact. ``shared`` provides only
    a before-contact physical opportunity: it serves the boundary from the
    preceding contact, including environmental contacts. There is a first pre
    wash but no final post wash. Role-boundary hand resets remain separate;
    they neither count as washing nor reset the worker's behavioral process.
    ``zone_separate`` instead groups each fixed bed and its environmental cell
    into a care zone, with one opportunity on physical entry and one on exit.
    Contacts inside a zone keep their transmission draws but add no hygiene.
    A horizon ending inside a zone has no artificial terminal exit wash. This
    entry/exit episode assumption is inspired by patient-zone models; it does
    not represent all clinical hygiene indications or a validated site workflow.
    Inactive returning staff complete a pending zone exit, but passing another
    zone without active care does not generate an entry or transmission contact.
    These are explicit scenario laws, not empirically fitted hygiene models.
    """

    def __init__(self, num_HCWs, num_Patients, num_Goo, prob_patient_sick,
                 prob_new_patient, prob_transmission, isolation_factor,
                 cleaningDay, isolate_sick, icu_hcw_wash_rate,
                 outside_hcw_wash_rate, height, width, *, seed=None,
                 stay_extension_days=0, los_mode='synthetic',
                 susceptible_los_mean_days=5., colonized_los_mean_days=14.,
                 data_collection_interval=None,
                 initial_environment_contaminated=True,
                 hygiene_opportunities='separate',
                 hygiene_autocorrelation=0., hygiene_efficacy=1.):
        super().__init__(seed=seed)
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        if any(not isinstance(value, numbers.Integral) or isinstance(value, bool)
               for value in (num_HCWs, num_Patients, num_Goo, width, height)):
            raise ValueError('Fixed topology counts and dimensions must be integers')
        if (num_HCWs, num_Patients, num_Goo, width, height) != (10, 30, 30, 32, 11):
            raise ValueError('Supported fixed topology requires 10 nurses, '
                             '30 patients/beds, 30 environment sites, width=32, height=11')
        for key, value in dict(prob_patient_sick=prob_patient_sick,
                               prob_new_patient=prob_new_patient,
                               prob_transmission=prob_transmission,
                               isolation_factor=isolation_factor,
                               icu_hcw_wash_rate=icu_hcw_wash_rate,
                               outside_hcw_wash_rate=outside_hcw_wash_rate).items():
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f'{key} must be finite and in [0, 1]')
        if not math.isfinite(cleaningDay) or cleaningDay <= 0:
            raise ValueError('cleaningDay must be positive and finite')
        if not math.isfinite(stay_extension_days) or stay_extension_days < 0:
            raise ValueError('stay_extension_days must be nonnegative and finite')
        if los_mode not in ('synthetic', 'historical'):
            raise ValueError("los_mode must be 'synthetic' or 'historical'")
        if hygiene_opportunities not in ('separate', 'shared', 'zone_separate'):
            raise ValueError("hygiene_opportunities must be 'separate', 'shared' or 'zone_separate'")
        if (not isinstance(hygiene_autocorrelation, numbers.Real)
                or not math.isfinite(hygiene_autocorrelation)
                or not 0 <= hygiene_autocorrelation < 1):
            raise ValueError('hygiene_autocorrelation must be finite and in [0, 1)')
        if (not isinstance(hygiene_efficacy, numbers.Real)
                or not math.isfinite(hygiene_efficacy)
                or not 0 <= hygiene_efficacy <= 1):
            raise ValueError('hygiene_efficacy must be finite and in [0, 1]')
        for value in (susceptible_los_mean_days, colonized_los_mean_days):
            if not math.isfinite(value) or value < 1:
                raise ValueError('geometric LOS means must be finite and at least one day')

        self.tick = 0
        self.total_ticks = 0
        self.ticks_in_hour = 108
        self.ticks_in_day = 24 * self.ticks_in_hour
        self.data_collection_interval = (self.ticks_in_day if data_collection_interval is None
                                         else data_collection_interval)
        if (not isinstance(self.data_collection_interval, numbers.Integral)
                or self.data_collection_interval < 1):
            raise ValueError('data_collection_interval must be a positive integer')
        self.stay_extension_days = stay_extension_days
        self.stay_extension_ticks = round(stay_extension_days * self.ticks_in_day)
        self.los_mode = los_mode
        self.susceptible_los_mean_days = susceptible_los_mean_days
        self.colonized_los_mean_days = colonized_los_mean_days
        self.num_HCWs = num_HCWs
        self.num_Patients = num_Patients
        self.num_Goo = num_Goo
        self.schedule = BaseScheduler(self)
        self.grid = MultiGrid(width, height, torus=False)
        self.prob_patient_sick = prob_patient_sick
        self.prob_new_patient = prob_new_patient
        self.prob_transmission = prob_transmission
        self.isolation_factor = isolation_factor
        self.cleaningDay = cleaningDay
        self.isolate_sick = isolate_sick
        self.icu_hcw_wash_rate = icu_hcw_wash_rate
        self.outside_hcw_wash_rate = outside_hcw_wash_rate
        self.hygiene_opportunities = hygiene_opportunities
        self.hygiene_autocorrelation = float(hygiene_autocorrelation)
        self.hygiene_efficacy = float(hygiene_efficacy)
        self._legacy_hygiene = (hygiene_opportunities == 'separate'
                               and hygiene_autocorrelation == 0
                               and hygiene_efficacy == 1)
        # Effective clears include successful clearing while already clean;
        # decontaminations count only actual dirty-to-clean hand transitions.
        self.hygiene_counter_names = ('hygiene_opportunity_count',
                                     'hygiene_completion_count',
                                     'hygiene_effective_clear_count',
                                     'hygiene_decontamination_count')
        for counter in self.hygiene_counter_names:
            setattr(self, counter, 0)
        self.hygiene_counts_by_role = {}
        self.initial_environment_contaminated = bool(initial_environment_contaminated)
        self.discharged = []
        self.current_patients = []
        self.empty_beds = set()
        self.shared_beds = []
        self.isol_beds = []
        # Zones belong to fixed bed locations, never patient IDs. A relocated
        # patient automatically belongs to the zone of its current bed.
        self.care_zone_by_position = {}
        self.nurse_list = []
        self.summon = False
        self.summoner = 0
        self.cumul_patients = 0
        self.cumul_sick_patients = 0
        self.num_infecByHCW = 0
        self.acquired_colonizations = 0
        self.imported_colonizations = 0
        self.discharged_patients = 0
        self.bed_visits = 0
        self.patient_contacts = 0
        self.environment_contacts = 0
        self.patient_days = 0.
        self.susceptible_patient_days = 0.
        self._patient_ticks = 0
        self._susceptible_patient_ticks = 0

        # Fixed activation order is retained and explicitly documented.
        for uid, x, y, care, shifts in ((-17, 31, 0, 30, 2),
                                      (-18, 30, 0, 20, 1),
                                      (-19, 30, 1, 3, 1)):
            self._place(XrayDr(uid, self, False, outside_hcw_wash_rate,
                               x, y, workHours=24, numCare=care, shiftsPerDay=shifts))
        for uid, x, y, hours in ((-16, 31, 10, 24), (-15, 30, 10, 12),
                                 (-14, 30, 9, 12), (-13, 31, 9, 3),
                                 (-12, 28, 10, 1), (-11, 28, 9, 1)):
            self._place(Dr(uid, self, False, icu_hcw_wash_rate,
                           x, y, workHours=hours, numCare=30))
        for uid in range(-num_HCWs, 0):
            nurse = Nurse(uid, self, False, icu_hcw_wash_rate, -2*uid-2, 0)
            nurse.path = paths[-uid-1]
            nurse.hall = 6 if nurse.path[0][0] == 'A' else 5
            self.nurse_list.append(nurse)
            self._place(nurse)

        for i in range(30):
            if i < 11:
                x, y = 2*i+5, 1
            elif i < 13:
                x, y = 2*(i-11)+1, 3
            elif i < 15:
                x, y = 2*(i-13)+27, 3
            elif i < 17:
                x, y = 2*(i-15)+1, 8
            elif i < 19:
                x, y = 2*(i-17)+27, 8
            else:
                x, y = 2*(i-19)+5, 10
            isolated = x in {1, 3, 5, 7, 23, 25, 27, 29}
            bed = (IsolatedBed if isolated else Bed)(i, self, False, x, y)
            self._place(bed)
            (self.isol_beds if isolated else self.shared_beds).append(bed)
            self.admit_patient(bed)
            environment_position = (x, y-1 if y > 5 else y+1)
            self.care_zone_by_position[(x, y)] = bed.unique_id
            self.care_zone_by_position[environment_position] = bed.unique_id
            self._place(Goo(-30-i, self, self.initial_environment_contaminated,
                            *environment_position))

        self.refresh_beds()
        self.datacollector = DataCollector(model_reporters={
            'Total ticks': 'total_ticks',
            'Number of Patients sick': getNumSick,
            'Number of HCW colonized': getNumColonized,
            'Cumulative number of colonized patients': getCumulNumSick,
            'HCW related infections': getHCWInfec,
            'Cumulative Patients': getCumul,
            'Imported colonizations': 'imported_colonizations',
            'Acquired colonizations': 'acquired_colonizations',
            'Patient days': 'patient_days',
            'Susceptible patient days': 'susceptible_patient_days',
            'Bed visits': 'bed_visits',
            'Patient contacts': 'patient_contacts',
            'Environment contacts': 'environment_contacts',
        })
        self.running = True
        self.datacollector.collect(self)
        self.assert_invariants()

    def _place(self, agent):
        self.schedule.add(agent)
        self.grid.place_agent(agent, (agent.x, agent.y))

    def make_hygiene_rng(self, worker_id):
        """Independent worker stream; creating it consumes no model RNG draws.

        An integer model seed and worker ID replay the same stream. Unseeded
        models use fresh entropy. This separation is not a cross-arm coupling:
        opportunities and epidemiological event draws may differ between arms.
        """
        entropy = (None if self.seed is None
                   else [int(self.seed), int(worker_id) & 0xffffffff, 0x485947])
        return np.random.default_rng(np.random.SeedSequence(entropy))

    def sample_stay_ticks(self, colonized):
        if colonized and self.los_mode == 'historical':
            # An unverified historical aggregate array, not a certified dataset.
            from routes import length_stay
            days = self.rng.choice(length_stay)
        else:
            mean = (self.colonized_los_mean_days if colonized
                    else self.susceptible_los_mean_days)
            days = self.rng.geometric(1. / mean)
        return int(days) * self.ticks_in_day

    def admit_patient(self, bed):
        patient = Patient(self.num_Patients + self.cumul_patients, self,
                          False, bed.x, bed.y)
        self._place(patient)
        self.cumul_patients += 1
        if patient.colonized:
            self.imported_colonizations += 1
            self.cumul_sick_patients += 1
        return patient

    def get_nurse_list(self):
        return list(self.nurse_list)

    def get_patient_list(self):
        return list(self.current_patients)

    def refresh_beds(self):
        self.empty_beds.clear()
        for bed in self.shared_beds + self.isol_beds:
            bed.checkFilled()
            if not bed.filled:
                self.empty_beds.add(bed.pos)

    def _remove_discharged(self):
        for patient in self.discharged:
            self.grid.remove_agent(patient)
            self.schedule.remove(patient)
            self.current_patients.remove(patient)
            self.discharged_patients += 1
        self.discharged.clear()

    def _move_patient(self, patient, destination):
        self.grid.move_agent(patient, destination)
        patient.x, patient.y = destination
        patient.checkIsolated()

    def _relocate_detected_patients(self):
        for bed in self.shared_beds:
            occupants = [a for a in self.grid.get_cell_list_contents([bed.pos])
                         if a.isPatient and a.positive]
            if not occupants:
                continue
            sick = occupants[0]
            for isolated_bed in self.isol_beds:
                isolated_bed.checkFilled()
                if isolated_bed.filledSick:
                    continue
                other = next((a for a in self.grid.get_cell_list_contents([isolated_bed.pos])
                              if a.isPatient), None)
                self._move_patient(sick, isolated_bed.pos)
                if other is not None:
                    self._move_patient(other, bed.pos)
                bed.checkFilled()
                isolated_bed.checkFilled()
                break

    def step(self):
        self._patient_ticks += len(self.current_patients)
        self._susceptible_patient_ticks += sum(not a.colonized for a in self.current_patients)
        self.patient_days = self._patient_ticks / self.ticks_in_day
        self.susceptible_patient_days = self._susceptible_patient_ticks / self.ticks_in_day
        self.total_ticks += 1
        self.tick = self.total_ticks % self.ticks_in_day
        if self.tick % self.ticks_in_hour == 0:
            self.summoner = self.tick // self.ticks_in_hour
            self.summon = 0 < self.summoner < 16
        self.schedule.step()
        self._remove_discharged()
        self.refresh_beds()
        if self.isolate_sick:
            self._relocate_detected_patients()
            self.refresh_beds()
        if self.total_ticks % self.data_collection_interval == 0:
            self.datacollector.collect(self)

    def assert_invariants(self):
        """Check patient conservation and agreement of all location records."""
        scheduled = [a for a in self.schedule.agents if a.isPatient]
        assert set(scheduled) == set(self.current_patients)
        assert len(scheduled) == len(set(scheduled))
        assert self.cumul_patients - self.discharged_patients == len(scheduled)
        assert self.cumul_sick_patients == self.imported_colonizations + self.acquired_colonizations
        assert self.num_infecByHCW == self.acquired_colonizations
        grid_patients = [a for cell, _, _ in self.grid.coord_iter()
                         for a in cell if a.isPatient]
        assert set(grid_patients) == set(scheduled)
        assert len(grid_patients) == len(scheduled)
        bed_positions = {bed.pos for bed in self.shared_beds + self.isol_beds}
        occupied = set()
        for patient in scheduled:
            assert patient.pos in bed_positions and patient.pos not in occupied
            occupied.add(patient.pos)
            assert (patient.x, patient.y) == patient.pos
            assert patient.isolated == (patient.pos[0] <= 7 or patient.pos[0] >= 23)
            assert patient.stay > 0
            assert patient in self.grid.get_cell_list_contents([patient.pos])
        assert self.empty_beds == bed_positions - occupied
        for bed in self.shared_beds + self.isol_beds:
            occupants = [a for a in self.grid.get_cell_list_contents([bed.pos]) if a.isPatient]
            assert bed.filled == bool(occupants)
            assert bed.filledSick == any(a.positive for a in occupants)

    def run_model(self, n):
        if not isinstance(n, numbers.Integral) or n < 0:
            raise ValueError('n must be a nonnegative integer tick count')
        for _ in range(n):
            self.step()
