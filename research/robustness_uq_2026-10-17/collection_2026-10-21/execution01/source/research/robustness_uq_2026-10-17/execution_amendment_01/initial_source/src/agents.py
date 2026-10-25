from mesa import Agent


ICUA = ['A14','A15','A1','A2','A3','A4','A5','A6','A7','A8','A9','A10','A11','A12','A13']
ICUB = ['B14','B15','B1','B2','B3','B4','B5','B6','B7','B8','B9','B10','B11','B12','B13']

###
def room2coord(room): # Outputs
        
        if len(room) == 2:
            pos = int(room[1]) # number after the letter
        else: # two digit number
            pos = int(room[1])*10 + int(room[2])
        
        letter = room[0] # letter
        if letter == 'A':
            if pos < 12: 
                y = 10
                x = 2 * abs(11-pos) + 5
            elif pos < 14:
                y = 8
                x = 2 * abs(pos-12) + 1
            else:
                y = 8
                x = 2 * abs(pos-14) + 27
        else: # 'B'
            if pos < 12: 
                y = 1
                x = 2 * abs(11-pos) + 5
            elif pos < 14:
                y = 3
                x = 2 * abs(pos-12) + 1
            else:
                y = 3
                x = 2 * abs(pos-14) + 27
        return x,y

class CPE_Agent(Agent):
    """ An agent with fixed colonized status."""
    def __init__(self, unique_id, model, colonized, x, y):
        super().__init__(unique_id, model)
        self.colonized = colonized
        self.isCPE_Agent = True
        self.isPatient = False
        self.isHCW = False
        self.isEnvironment = False
        self.isBed = False
        self.isGoo = False
        self.isWall = False
        self.x = x
        self.y = y

class HCW(CPE_Agent):
    """ An agent with fixed colonized status. Can wash hands"""
    def __init__(self, unique_id, model, colonized, hand_wash_rate, x, y): #path):
        super().__init__(unique_id, model, colonized, x, y)
        self.hand_wash_rate = hand_wash_rate
        self._hygiene_previous_completion = None
        self._hygiene_zone = None
        self.hygiene_rng = (None if model._legacy_hygiene
                            else model.make_hygiene_rng(unique_id))
        for counter in model.hygiene_counter_names:
            setattr(self, counter, 0)
        self.hygiene_role = type(self).__name__
        model.hygiene_counts_by_role.setdefault(
            self.hygiene_role, dict.fromkeys(model.hygiene_counter_names, 0))
        
        self.test_rate = 1
        self.true_positive_rate = 1
        
        self.isHCW = True
        self.moveTick = 0
        self.path = []
        self.reached = False
        self.reachedHall = False
        self.hall = 5
        self.bed_visits = 0
        self.patient_contacts = 0
        self.environment_contacts = 0
        self.last_contact_tick = None
        self.last_contact_pos = None
        self.current_destination = None
        self.activity = 'ready'
        self._journey = None
        self._arrival_recorded = False
        self._route_entry_completed = False
        self._contact_position = None
        self._contacted_objects = set()

    def _reset_journey(self):
        self._journey = None
        self.reachedHall = False
        self._arrival_recorded = False

    def _begin_care_episode(self):
        # A repeated route entry is a new task on a subsequent tick. Merely
        # waiting at the same patient/site is still the same contact episode.
        self._contacted_objects.clear()

    def _travel(self, destination, activity):
        """Take at most one orthogonal cell via the assigned hallway.

        The recovered layout supplies no wall/door geometry. If a Wall agent
        is supplied, a blocked step waits; it never walks through that wall.
        """
        destination = tuple(destination)
        journey = (destination, activity)
        if journey != self._journey:
            same_arrival = (self._journey is not None
                            and self._journey[0] == destination
                            and self.pos == destination)
            self._journey = journey
            self.reachedHall = self.pos[1] == self.hall
            if not same_arrival:
                self._arrival_recorded = False
        self.current_destination = destination
        self.activity = activity

        if self.pos != destination:
            x, y = self.pos
            if not self.reachedHall:
                new_position = (x, y + (1 if y < self.hall else -1))
            elif x != destination[0]:
                new_position = (x + (1 if x < destination[0] else -1), y)
            else:
                new_position = (x, y + (1 if y < destination[1] else -1))
            if (self.model.grid.out_of_bounds(new_position)
                    or any(a.isWall for a in
                           self.model.grid.get_cell_list_contents([new_position]))):
                self.activity = 'blocked'
                return False
            self.model.grid.move_agent(self, new_position)
            self._contacted_objects.clear()
            # Check the updated coordinate: no extra hallway no-op tick.
            if self.pos[1] == self.hall:
                self.reachedHall = True

        arrived = self.pos == destination
        if arrived:
            self.activity = 'off duty' if activity == 'returning' else (
                'waiting for team' if activity == 'summoned' else 'at bed')
            if not self._arrival_recorded and activity != 'returning':
                if any(a.isBed for a in self.model.grid.get_cell_list_contents([self.pos])):
                    self.bed_visits += 1
                    self.model.bed_visits += 1
                self._arrival_recorded = True
        return arrived

    def move(self):
        """Visit one route entry, retaining its arrival for this tick's care."""
        if not len(self.path):
            self.current_destination = None
            self.activity = 'no route'
            return False
        self.moveTick %= len(self.path)
        if self._route_entry_completed:
            self.moveTick = (self.moveTick + 1) % len(self.path)
            self._route_entry_completed = False
            self._reset_journey()
            self._begin_care_episode()
        arrived = self._travel(room2coord(self.path[self.moveTick]), 'rounds')
        self._route_entry_completed = arrived
        return arrived

    def answerSummon(self, destination):
        return self._travel(destination, 'summoned')

        
    def testCRE(self, other):
        test = self.model.rng.choice([1,0], p = [self.test_rate, 1-self.test_rate])
        if test:
            positive = self.model.rng.choice([1,0], p = [self.true_positive_rate, 1-self.true_positive_rate])
            if positive:
                other.positive = True

    def spread(self):
        """Apply declared physical hygiene opportunities per care episode.

        This is idealized immediate screening; carriage and clinical infection
        are not interchangeable endpoints. A continuous stay does not repeat
        care, transmission draws or hygiene merely because another tick passed.
        """
        if self.model.hygiene_opportunities == 'zone_separate':
            self._update_zone_hygiene()
        self.reached = False
        eligible = [a for a in self.model.grid.get_cell_list_contents([self.pos])
                    if a.isPatient or a.isGoo]
        if self._contact_position != self.pos:
            self._contacted_objects.clear()
        self._contact_position = self.pos
        present = {a.unique_id for a in eligible}
        self._contacted_objects.intersection_update(present)
        for other in eligible:
            if other.unique_id in self._contacted_objects:
                continue
            self._contacted_objects.add(other.unique_id)
            self.reached = True
            counter = 'patient_contacts' if other.isPatient else 'environment_contacts'
            setattr(self, counter, getattr(self, counter) + 1)
            setattr(self.model, counter, getattr(self.model, counter) + 1)
            self.last_contact_tick = self.model.total_ticks
            self.last_contact_pos = self.pos
            if self.model.hygiene_opportunities != 'zone_separate':
                self.handwash()
            prob = self.model.prob_transmission
            if other.isolated:
                prob *= self.model.isolation_factor
            if self.colonized and not other.colonized:
                if self.model.rng.random() < prob:
                    other.colonized = True
                    if other.isPatient:
                        other.state = 'colonized'
                        other.stay += self.model.stay_extension_ticks
                        self.model.cumul_sick_patients += 1
                        self.model.num_infecByHCW += 1
                        self.model.acquired_colonizations += 1
            elif not self.colonized and other.colonized:
                if self.model.rng.random() < prob:
                    self.colonized = True
            # Detection is separate from the worker's contamination state.
            if other.isPatient and other.colonized:
                self.testCRE(other)
            if self.model.hygiene_opportunities == 'separate':
                self.handwash()

    def _update_zone_hygiene(self, *, allow_entry=True):
        """One physical pre/post opportunity on actual care-zone changes.

        The bed and its environment share a fixed zone. Waiting, route task
        advancement and activity changes inside it do not create new boundaries.
        Moving directly between zones executes old-zone exit then new-zone entry.
        Leaving for a hallway executes exit and records no current zone. Implicit
        role hand resets are separate events and do not reset behavioral memory.
        Inactive returning movement only completes a pending exit; merely
        passing through another zone off duty is not a new care episode.
        """
        zone = self.model.care_zone_by_position.get(self.pos)
        if not allow_entry:
            if self._hygiene_zone is not None and zone != self._hygiene_zone:
                self.handwash()
                self._hygiene_zone = None
            return
        if zone == self._hygiene_zone:
            return
        if self._hygiene_zone is not None:
            self.handwash()
        if zone is not None:
            self.handwash()
        self._hygiene_zone = zone

    def handwash(self):
        """Advance behavior even with clean hands, then attempt binary clearing.

        Completion is stationary at initialization. Its Markov transition is
        h + rho*(1-h) after completion, h*(1-rho) after failure. Shift resets
        clear hands independently and leave this behavioral state unchanged.
        The legacy configuration consumes exactly its original one RNG draw.
        Other configurations use only the dedicated stream for hygiene.
        """
        self._count_hygiene('hygiene_opportunity_count')
        if self.model._legacy_hygiene:
            completed = self.model.rng.random() < self.hand_wash_rate
        else:
            probability = self.hand_wash_rate
            rho = self.model.hygiene_autocorrelation
            if self._hygiene_previous_completion is not None:
                probability = (probability + rho * (1 - probability)
                               if self._hygiene_previous_completion
                               else probability * (1 - rho))
            completed = self.hygiene_rng.random() < probability
        self._hygiene_previous_completion = bool(completed)
        if not completed:
            return
        self._count_hygiene('hygiene_completion_count')
        efficacy = self.model.hygiene_efficacy
        effective = (efficacy == 1 or
                     (efficacy > 0 and self.hygiene_rng.random() < efficacy))
        if effective:
            self._count_hygiene('hygiene_effective_clear_count')
            if self.colonized:
                self._count_hygiene('hygiene_decontamination_count')
            self.colonized = False

    def _count_hygiene(self, counter):
        setattr(self, counter, getattr(self, counter) + 1)
        setattr(self.model, counter, getattr(self.model, counter) + 1)
        self.model.hygiene_counts_by_role[self.hygiene_role][counter] += 1
            
class Nurse(HCW):
    def __init__(self, unique_id, model, colonized, hand_wash_rate, x, y): #path):
        super().__init__(unique_id, model, colonized, hand_wash_rate, x, y)
        self.workHours = 8

    def step(self):
        # A new eight-hour role slot starts clean before its first contact.
        if self.model.tick % (self.model.ticks_in_hour * self.workHours) == 0:
            self.colonized = False
        if self.model.summon:
            if self.model.summoner > 0 and self.model.summoner < 16: 
                if self.model.summoner < 8 or self.model.summoner > 13: # on the east wing of the ICU
                    nurses_A = {-1,-2,-3} # ID of the nurses in east wing
                    nurses_B = {-6,-7,-8}
                
                else:
                    nurses_A = {-3,-4,-5} # ID of the nurses in west wing
                    nurses_B = {-10,-9,-8}
                
                destination = str(self.model.summoner) # bed in need of position change
                summoner_loc_A = self.model.grid.get_cell_list_contents([room2coord('A' + destination)]) # location of bed in need of position change
                summoner_loc_B = self.model.grid.get_cell_list_contents([room2coord('B' + destination)]) # location of bed in need of position change
                # turning off model.summon when all 3 nurses are present
                ppl_at_loc_A = set()
                ppl_at_loc_B = set()
                for x in summoner_loc_A:
                    ppl_at_loc_A.add(x.unique_id)
                for x in summoner_loc_B:
                    ppl_at_loc_B.add(x.unique_id)


                # if you are called
                if self.unique_id in nurses_A:
                    self.answerSummon(room2coord('A' + destination)) # go to this room
                    self.spread()
                elif self.unique_id in nurses_B:
                    self.answerSummon(room2coord('B' + destination)) # go to this room
                    self.spread()
                # if not called
                else:
                    self.move()
                    self.spread()
                
                if nurses_A.issubset(ppl_at_loc_A) and nurses_B.issubset(ppl_at_loc_B): # all 3 nurses have joined
                    for x in summoner_loc_A + summoner_loc_B:
                        if isinstance(x, Nurse):
                            x.reachedHall = False
                    self.model.summon = False
            
        # regular activity
        else:
            self.move()
            self.spread()
            
        
        

            
            
class Dr(HCW):
    def __init__(self, unique_id, model, colonized, hand_wash_rate, x, y, workHours, numCare): #path):
        super().__init__(unique_id, model, colonized, hand_wash_rate, x, y)
        self.original_pos = (x,y)
        self.workHours = workHours
        self.activated = False
        self.path = self.model.random.sample(ICUA+ICUB, numCare)
        self.startTime = 0
        self.endTime = self.startTime + self.model.ticks_in_hour * self.workHours

    def leaveICU(self):
        arrived = self._travel(self.original_pos, 'returning')
        if self.model.hygiene_opportunities == 'zone_separate':
            self._update_zone_hygiene(allow_entry=False)
        return arrived
    
    def step(self):
        # Role slots restart daily; 24-hour slots are active across midnight.
        phase = self.model.tick
        if phase in (self.startTime, self.endTime % self.model.ticks_in_day):
            self.colonized = False
            self._reset_journey()
        self.activated = self.startTime <= phase < self.endTime
        if self.activated:
            self.move()
            self.spread()
        else:
            self.leaveICU()

class XrayDr(Dr):
    def __init__(self, unique_id, model, colonized, hand_wash_rate, x, y, workHours, numCare, shiftsPerDay):
        super().__init__(unique_id, model, colonized, hand_wash_rate, x, y, workHours, numCare)
        self.shiftsPerDay = shiftsPerDay    
        self.original_pos = (x,y)
        # The model schedules its first completed tick as 1, so a phase-zero
        # round must start active rather than wait until the next day boundary.
        self.activated = True
        self.path = self.model.random.sample(ICUA+ICUB, numCare)
        self.k = len(self.path)
        self.last_patient = self.path[numCare-1]
    
    
    def step(self):
        if self.model.tick % (self.model.ticks_in_day / self.shiftsPerDay)== 0:
            self.colonized = False
            self.activated = True
            self.moveTick = 0
            self._route_entry_completed = False
            self._reset_journey()
        # Finish the final route entry only after its arrival-tick contact.
        # Coordinate coincidence with last_patient is not round completion.
        if (self.activated and self._route_entry_completed
                and self.moveTick == len(self.path) - 1):
            self.activated = False
            self.colonized = False
            self.moveTick = 0
            self._route_entry_completed = False
            self._reset_journey()
        
        if self.activated:
            self.move()
            self.spread()
        else:
            self.leaveICU()
        
class Patient(CPE_Agent):
    """ An agent with colonized status.
    Cannot move and leaves hospital after a certain time"""
    def __init__(self, unique_id, model, colonized, x, y):
        super().__init__(unique_id, model, colonized, x, y)
        self.isPatient = True
        #self.infecByHCW = False # in the beginning, nobody is infected by HCW
        
        self.checkIsolated()
        self.model.current_patients.append(self)
        #self.state
        sick = self.model.rng.choice([1,2], p=[self.model.prob_patient_sick,1-self.model.prob_patient_sick])
        if sick == 1:
            self.colonized = True
        self.state = 'colonized' if self.colonized else 'susceptible'

        
        self.stay = self.model.sample_stay_ticks(self.colonized)

        #self.new = True # SCV ready to go sir

        self.positive = False
        
    def checkIsolated(self):
        xpos = self.pos[0] if self.pos is not None else self.x
        if xpos <= 7 or xpos >= 23:
            self.isolated = True
        else:
            self.isolated = False

    def step(self):
        #remove oneself if the stay is too long
        self.stay -= 1
        if self.stay <= 0 and self not in self.model.discharged:
            self.model.discharged.append(self)
        
        #if self.new == True:
        #    print("New Patient {}: {}, {} days left".format(self.unique_id, self.state, self.stay))
        #else:
        #    print("Patient {}: {}, {} days left".format(self.unique_id, self.state, self.stay))
        #self.new = False

        self.checkIsolated() #alters state if moved to/from isolated bed
    
    
#%%


class Environment(CPE_Agent):
    """ An agent with colonized status.
    Cannot move and does not leave."""
    def __init__(self, unique_id, model, colonized, x, y):
        super().__init__(unique_id, model, colonized, x, y)
        self.isEnvironment = True    
        self.filled = False # all patients fill the beds
        self.filledSick = False
        self.isGoo = False
        self.isBed = False
        self.isIsolatedBed = False
        self.isWall = False
        



class Wall(Environment):
    """ An agent with colonized status.
    HCW cannot pass through walls."""
    def __init__(self, unique_id, model, colonized, x, y):
        super().__init__(unique_id, model, colonized, x, y)
        self.isWall = True
    def step(self):
        pass

class Bed(Environment):
    """ An agent with colonized status.
    Can get contaminated, but cleans after a patient leaves."""
    def __init__(self, unique_id, model, colonized, x, y):
        super().__init__(unique_id, model, colonized, x, y)
        self.isBed = True
        self.isIsolatedBed = False
    
    def step(self):
        self.checkFilled()
        if not self.filled and self.model.rng.random() < self.model.prob_new_patient:
            self.model.admit_patient(self)
        self.checkFilled()
        if self.filled:
            self.model.empty_beds.discard(self.pos)
        else:
            self.model.empty_beds.add(self.pos)

    def checkFilled(self):
        self.filled = False
        self.filledSick = False #reset
        cellmates = self.model.grid.get_cell_list_contents([self.pos])
        if len(cellmates)>1:
            for cellmate in cellmates:
                if cellmate.isPatient:
                    self.filled = True
                    if self.filled and cellmate.positive:
                        self.filledSick = True
                    break # other contents of bed doesn't matter
    
        
    
    
class IsolatedBed(Bed):
    """ An agent with colonized status.
    Can get contaminated, but cleans after a patient leaves.
    The isolated bed has lower probability of transmission compared to other beds."""
    def __init__(self, unique_id, model, colonized, x, y):
        super().__init__(unique_id, model, colonized, x, y)
        self.isIsolatedBed = True
        
class Goo(Environment):
    """ An agent with colonized status.
    Can get contaminated, cleans once in a while."""
    def __init__(self, unique_id, model, colonized, x, y):
        super().__init__(unique_id, model, colonized, x, y)
        self.isGoo = True
        self.clean_tick = self.model.cleaningDay * self.model.ticks_in_day
        
        if self.x <= 7 or self.x >= 23:
            self.isolated = True
        else:
            self.isolated = False
    def handwash(self): # 100%
        self.colonized = False
        
    def step(self):
        #self.checkFilled()
        self.clean_tick -= 1
        
        if (self.clean_tick <= 0):
            self.handwash()
            self.clean_tick = self.model.cleaningDay * self.model.ticks_in_day
        
    

# make a patient class, test run
