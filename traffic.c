#include types.h
#include lib.h
#include test.h
#include clock.h
#include thread.h
#include synch.h
#include synchprobs.h
#include lamebusltimer.h
#include kernerrno.h


#define MAX_THREADS 10

static int NumIterations = 100;   number of vehicle arrivals per thread
static int NumThreads = 10;       number of concurrent simulation threads
static int InterArrivalTime = 1;       time between vehicle arrivals
static int ServiceTime = 1;     time in the intersection
static int DirectionBias = 0;   0 = unbiased, 1 = biased

static struct semaphore SimulationWait;


typedef struct Vehicles
{
  Direction origin;
  Direction destination;
} Vehicle;

  vehicles in the intersection 
  note this is an array of volatile pointers 
static Vehicle  volatile vehicles[MAX_THREADS];
 semaphore to protect vehicles array 
static struct semaphore mutex;

 performance statistics, one set per origin direction 
static volatile time_t total_wait_secs[4];
static volatile uint32_t total_wait_nsecs[4];
static volatile int wait_count[4];
static volatile time_t max_wait_secs[4];
static volatile uint32_t max_wait_nsecs[4];
 mutex to provide mutual exclusion to performance stats 
static struct semaphore perf_mutex;
 number of milliseconds per timer tick (LT_GRANULARITY is in usec) 
#define TICK_MILLISEC (LT_GRANULARITY1000)

 simulation start and end time 
time_t start_sec, end_sec;
uint32_t start_nsec, end_nsec;

 bias direction, for arrival biasing 
Direction heavy_direction;

 functions defined and used internally 
static void initialize_state(void);
static void cleanup_state(void);
static void print_direction(Direction d);
static void print_perf_stats(void);
static void vehicle_simulation(void ptr, unsigned long thread_num);
static bool right_turn(Vehicle v);
static void check_constraints(int thread_num);

static void
print_direction(Direction d) {
  switch (d)
    {
    case north
      kprintf(N);
      break;
    case east
      kprintf(E);
      break;
    case south
      kprintf(S);
      break;
    case west
      kprintf(W);
      break;
    }
}    

static Direction
choose_direction(void) {
  int x;
  if (DirectionBias) {
    x = random()%10;
    if (x = 4) {
      return heavy_direction;
    }
    else {
      return x;
    }
  }
  else {
    return random()%4;
  }
}    

static void
print_perf_stats(void) {
  int i;
  int max_allowed_ms;
  int wait_msecs,mean_wait_msecs,max_wait_msecs;
  int total_wait_msecs = 0;
  int total_count = 0;
  int sim_msec;
  time_t run_sec;
  uint32_t run_nsec;

  max_allowed_ms = (NumThreads-1)ServiceTimeTICK_MILLISEC;
  kprintf(max permitted waitt%d.%03d secondsn,
	  max_allowed_ms1000,max_allowed_ms%1000);
  for(i=0;i4;i++) {
    print_direction((Direction)i);
    kprintf(t);
    if (wait_count[i]  0) {
      wait_msecs = (total_wait_secs[i]1000+total_wait_nsecs[i]1000000);
      total_wait_msecs += wait_msecs;
       some rounding error here, in millisecond range
      mean_wait_msecs = wait_msecswait_count[i];
      total_count += wait_count[i];
      max_wait_msecs = max_wait_secs[i]1000+max_wait_nsecs[i]1000000;
      kprintf(%d vehicles, average wait %d.%03d seconds, max wait %d.%03d secondsn,
	      wait_count[i], mean_wait_msecs1000,mean_wait_msecs%1000,
	      max_wait_msecs1000,max_wait_msecs%1000);
    } else {
      kprintf(0 vehicles, average wait 0.000 seconds, max wait 0.000 secondsn);
    }
  }
   then the average wait time for all vehicles 
  if (total_count  0) {
    kprintf(allt%d vehicles, average %d.%03d seconds waitingn,total_count,
	    (total_wait_msecstotal_count)1000,
	    (total_wait_msecstotal_count)%1000);
  } else{
    kprintf(allt0 vehicles, average 0.000 seconds waitingn);
  }
   finally, overall simulation run-time and throughput 
  getinterval(start_sec,start_nsec,end_sec,end_nsec,&run_sec,&run_nsec);
  sim_msec = run_sec1000;
  sim_msec += run_nsec1000000;
  kprintf(Simulation %d.%03d seconds, %d vehiclesn,
	  sim_msec1000,
	  sim_msec%1000,
	  total_count);
} 

bool
right_turn(Vehicle v) {
  KASSERT(v != NULL);
  if (((v-origin == west) && (v-destination == south)) 
      ((v-origin == south) && (v-destination == east)) 
      ((v-origin == east) && (v-destination == north)) 
      ((v-origin == north) && (v-destination == west))) {
    return true;
  } else {
    return false;
  }
}

void
check_constraints(int thread_num) {
  int i;
  KASSERT(thread_num  NumThreads);
   compare newly-added vehicle to each other vehicles in in the intersection 
  for(i=0;iNumThreads;i++) {
    if ((i==thread_num)  (vehicles[i] == NULL)) continue;
     no conflict if both vehicles have the same origin 
    if (vehicles[i]-origin == vehicles[thread_num]-origin) continue;
     no conflict if vehicles go in opposite directions 
    if ((vehicles[i]-origin == vehicles[thread_num]-destination) &&
        (vehicles[i]-destination == vehicles[thread_num]-origin)) continue;
     no conflict if one makes a right turn and 
       the other has a different destination 
    if ((right_turn(vehicles[i])  right_turn(vehicles[thread_num])) &&
	(vehicles[thread_num]-destination != vehicles[i]-destination)) continue;
    kprintf(Vehicle A );
    print_direction(vehicles[i]-origin);
    kprintf(-);
    print_direction(vehicles[i]-destination);
    kprintf(n);
    kprintf(Vehicle B );
    print_direction(vehicles[thread_num]-origin);
    kprintf(-);
    print_direction(vehicles[thread_num]-destination);
    kprintf(n);
    panic(intersection synchronization constraint violation!n);
  }
}


static void
initialize_state(void)
{
  int i;
  for(i=0;iMAX_THREADS;i++) {    
    vehicles[i] = (Vehicle  volatile)NULL;
  }
  for(i=0;i4;i++) {
    total_wait_secs[i] = total_wait_nsecs[i] = wait_count[i] = 0;
    max_wait_secs[i] = max_wait_nsecs[i] = 0;
  }
  mutex = sem_create(Vehicle Mutex,1);
  if (mutex == NULL) {
    panic(could not create vehicle mutex semaphoren);
  }
  perf_mutex = sem_create(PerfMutex,1);
  if (perf_mutex == NULL) {
    panic(could not create perf_mutex semaphoren);
  }
  SimulationWait = sem_create(SimulationWait,0);
  if (SimulationWait == NULL) {
    panic(could not create SimulationWait semaphoren);
  }
  heavy_direction = random()%4;
   initialization for synchronization code 
  intersection_sync_init();

}



static void
cleanup_state(void)
{
  sem_destroy(mutex);
  sem_destroy(perf_mutex);
  sem_destroy(SimulationWait);
  intersection_sync_cleanup();
}

static void
in_intersection(void) {
   clocknap(ServiceTime);
}



static
void
vehicle_simulation(void  unusedpointer, 
               unsigned long thread_num)
{
  int i;
  Vehicle v;
  time_t before_sec, after_sec, wait_sec;
  uint32_t before_nsec, after_nsec, wait_nsec;
  int sleeptime;
  (void) unusedpointer;

  KASSERT((long)thread_num  NumThreads);
  for(i=0;iNumIterations;i++) {

    sleeptime = InterArrivalTime + random()%3 - 1;
    KASSERT(sleeptime = InterArrivalTime-1);
    KASSERT(sleeptime = InterArrivalTime+1);
    clocknap(sleeptime);
    if (random()%NumThreads   thread_num) {
      thread_yield();
    }
    
     choose where this vehicle is coming from 
    v.origin = choose_direction();
     choose where this vehicle is heading 
    v.destination = v.origin + (random()%3) + 1;
    if (v.destination = 4) {
      v.destination = v.destination % 4;
    }
    KASSERT(4  v.origin);
    KASSERT(4  v.destination);
    KASSERT(v.origin != v.destination);

    gettime(&before_sec,&before_nsec);
    intersection_before_entry(v.origin, v.destination);
    gettime(&after_sec,&after_nsec);

    P(mutex);
    KASSERT(vehicles[thread_num] == NULL);
    vehicles[thread_num] = &v;

    check_constraints(thread_num);
    V(mutex);

    in_intersection();

    P(mutex);
    KASSERT(vehicles[thread_num] == &v);
    vehicles[thread_num] = NULL;
    V(mutex);

    intersection_after_exit(v.origin, v.destination);

    getinterval(before_sec,before_nsec,after_sec,after_nsec,&wait_sec,&wait_nsec);
    P(perf_mutex);
    total_wait_secs[v.origin] += wait_sec;
    total_wait_nsecs[v.origin] += wait_nsec;
    if (total_wait_nsecs[v.origin]  1000000000) {
      total_wait_nsecs[v.origin] -= 1000000000;
      total_wait_secs[v.origin] ++;
    }
    wait_count[v.origin]++;
    if (wait_sec  max_wait_secs[v.origin]) {
      max_wait_secs[v.origin] = wait_sec;
      max_wait_nsecs[v.origin] = wait_nsec;
    }
    else if ((wait_sec == max_wait_secs[v.origin])&&
	     (wait_nsec  max_wait_nsecs[v.origin])) {
      max_wait_nsecs[v.origin] = wait_nsec;
    }
    V(perf_mutex);
  }

   indicate that this simulation is finished 
  V(SimulationWait); 
}


int
traffic_simulation(int nargs,
		  char  args)
{
  int i;
  int error;

  if ((nargs != 1) && (nargs != 6)) {
    kprintf(Usage command [threads iterations interarrivaltime servicetime directionbiasn);
    return EINVAL;   return failure indication
  }

  if (nargs == 6) {
    NumThreads = atoi(args[1]);
    if (NumThreads = 0  NumThreads  MAX_THREADS) {
      kprintf(invalid number of threads %dn,NumThreads);
      return EINVAL;
    }
    NumIterations = atoi(args[2]);
    if (NumIterations  0) {
      kprintf(invalid number of iterations per thread %dn,NumIterations);
      return EINVAL;
    }
    InterArrivalTime = atoi(args[3]);
    if (InterArrivalTime  0) {
      kprintf(invalid interarrival time %dn,InterArrivalTime);
      return EINVAL;
    }
    ServiceTime = atoi(args[4]);
    if (ServiceTime  0) {
      kprintf(invalid service time %dn,ServiceTime);
      return EINVAL;
    }
    DirectionBias = atoi(args[5]);
    if ((DirectionBias != 0) && (DirectionBias != 1)) {
      kprintf(Direction Bias (%d) must be 0 (unbiased) or 1 (biased)n,DirectionBias);
      return EINVAL;
    }
  }
  
  kprintf(Threads %d Iterations %d Interarrival time %d  Service time %d  Bias %sn,
          NumThreads,NumIterations,InterArrivalTime,ServiceTime,DirectionBiasyesno);

   initialize our simulation state 
  initialize_state();

   get simulation start time 
  gettime(&start_sec,&start_nsec);

  for (i = 0; i  NumThreads; i++) {
    error = thread_fork(vehicle_simulation thread, NULL, vehicle_simulation, NULL, i);
    if (error) {
      panic(traffic_simulation thread_fork failed %sn, strerror(error));
    }
  }
  
   wait for all of the vehicle simulations to finish before terminating   
  for(i=0;iNumThreads;i++) {
    P(SimulationWait);
  }

   get simulation end time 
  gettime(&end_sec,&end_nsec);

   clean up the simulation state 
  cleanup_state();

   display performance stats 
  print_perf_stats();

  return 0;
}

