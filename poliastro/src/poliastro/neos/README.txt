2014-Aug-18                          DASTCOM5 
                                 USER QUICK START 

 CONTENTS

  BACKGROUND

  DISTRIBUTION RETRIEVAL

  USAGE
   Building software and library
   Running DXLOOK
   Using the FORTRAN reader library
    Example 1: Open DASTCOM5, return all information in a record
    Example 2: Open DASTCOM5, return only DASTCOM3 information
    Example 3: Open DASTCOM5, return only DASTCOM3 data in legacy D3READ order
    Example 4: Open legacy DASTCOM3, return data in legacy D3READ order

  PROGRAMMER NOTES
   Overview of database structure
    DAST5
    DCOM5
    ADDRESSING
    SCANNING THE DATABASE
    INDEX FILE 
   Alternative access methods

  CONTENTS OF DASTCOM5

  DASTCOM5 BYTE MAP
  
  CONTACT
  

BACKGROUND

 DASTCOM5 is a direct-access binary database. It contains heliocentric 
 ecliptic osculating elements for the known asteroids and comets, at single 
 instants, determined by a least-squares orbit solution fit to optical and 
 radar astrometric measurements. 

 Additional physical, dynamical, and covariance parameters are included when
 known. A total of 142 parameters per object are defined within DASTCOM5
 though many may not yet be assigned values, being undetermined. See the 
 "./dastcom5/src/curdef.inc" file for detailed information on the defined 
 fields, after retrieving and unzipping the distribution file described below. 

 This information is suitable for initializing numerical integrators, assessing
 orbit geometry, computing trajectory uncertainties, visual magnitude and, to 
 some extent, summarizing physical characterization of the body. It is the 
 database used by the JPL Horizons ephemeris system.

 Each asteroid has one data record containing the most recent information for 
 that object. Comets, however, may have multiple records, reflecting solutions 
 at different apparitions. This is because acceleration due to non-gravitational
 out-gassing can vary during different perihelion passages and be discernable 
 from the measurement data. Unless doing historical research, the record for 
 a comets' apparition closest to the present is typically most relevant.

 DASTCOM5 is a portable, subset product exported from a larger MySQL-based 
 relational master database called SBDB ("Small-Body Database") maintained at 
 the NASA/Jet Propulsion Laboratory.

 The DASTCOM5 distribution file archive is updated when necessary to include 
 newly discovered objects or orbit solution updates derived from newly reported
 measurements. Database updates typically occur a couple times per day, but 
 potentially as frequently as once per hour. They are placed on the public 
 FTP site of JPL's Solar System Dynamics Group (below).

 The DASTCOM5 distribution .zip file archive contains the following:

    * two binary database files (one holding asteroid data, the second holding
       comet data), 

    * a plain-text index (linking all objects to their DASTCOM5 record, 
       permitting look-up based on name, designation, SPK ID, packed MPC 
       designation, and historical aliases), 

    * documentation,

    * latest database reader source code (FORTRAN). The software has been 
       tested using GFORTRAN (gfortran), Lahey (lf95), Intel (ifort), and 
       SunStudio (f95) compilers in both 32 and 64-bit builds, under RedHat 
       Linux 4/5/6.

 The legacy DASTCOM3 file, in use over the last 20+ years (without an official 
 reader routine), will continue to be produced by conversion indefinately, 
 using DASTCOM5 as a source and discarding and reorganizing data into the 
 old structure. 

 However, DASTCOM3 cannot contain more than 999,999 object records and so 
 cannot be produced after that number of objects are cataloged. At the time of
 this writing, there are ~624,000 small-body records. While it is likely
 discoveries won't approach the DASTCOM3 limit for more than three years, it 
 is recommended users switch to DASTCOM5 well before that time. Note that 
 DASTCOM3 does not provide enough information to properly initialize modern 
 trajectory propagations for some objects having non-gravitational dynamics 
 (primarily comets with A3 parameters), and has therefore been obsolete for
 that purpose for several years.

 DASTCOM4 was an extension to DASTCOM3, defined in 2001, and primarily used 
 internally at JPL.

 The outline below assumes UNIX/Linux/Mac command line functionality, but 
 should be adaptable to PC/Windows environments.


DISTRIBUTION RETRIEVAL:

 1. Retrieve the latest database distribution by anonymous FTP from the
     JPL Solar System Dynamics server: 

     From command-line

       wget ftp://ssd.jpl.nasa.gov/pub/xfr/dastcom5.zip

     OR browser URL bar

       ftp://ssd.jpl.nasa.gov/pub/xfr/dastcom5.zip
    
     OR command-line FTP login

       ftp ssd.jpl.nasa.gov
       Name: anonymous
       Password: {your email address}
       cd pub/xfr
       binary
       get dastcom5.zip
       exit

    IMPORTANT NOTE: 
     The database could be updated at any moment within the interval 30-32 
     minutes after any hour. Make sure your retrieval is not occurring during 
     this time to avoid creating a local database on your system potentially 
     corrupted due to the remote source file being almost instantly rewritten 
     with an update while your slower Internet transfer is in progress.

 2. Unzip the archived directory structure 

       unzip -ao dastcom5.zip
 
     The "-ao" options instruct unzip to convert text files to local forms
     and over-write any existing files of the same name.

     A directory called "dastcom5" will be created in the current working
     directory, along with five sub-directories: 

       ./dastcom5/dat 
                     /dast5_le.dat  (binary asteroid data)
                     /dcom5_le.dat  (binary comet data)
                     /dastcom.idx   (plain-text database index)
                 /doc 
                     /README.txt
                     /news.txt
                     /dastcom5.map
                 /src               (FORTRAN source for libraries & apps)
                     /*.f 
                     /*.inc
                     /Makefile
                     /dxlook.inp
                 /exe               (initially empty)
                 /lib               (initially empty)
  
     The /dat sub-directory contains the latest matched-set of databases.

     The database files can be moved and renamed, but the rest of this 
     discussion assumes the directory structure above for the purpose of 
     example.

     Note that the three database files in a distribution make up a "matched 
     set" and MUST be replaced/updated together. If an old 'dast5' is used with
     a more recent 'dcom5' (or index), errors may result, including retrieval 
     of a different object than expected.

     Changes to the distribution or other developments will be noted in the 
     file './dastcom5/doc/news.txt'

USAGE

 It is recommended the reader library and applications in the './dastcom5/src'
 directory be used to read the database. Advantages include:

   - Binary-independent readability across IEEE platforms. It is not necessary
      to have a version of the DASTCOM database that is binary-compatible with 
      your computer; the reader software translates as necessary.

   - Legacy database compatibility; transparently supports prior DASTCOM3 and 
      DASTCOM4 databases, and the old limited-release D3READ reader subroutine,
      simplifying transition to the new database.

   - Transparent updates: if DASTCOM5 is altered in the future to include 
      additional information, user software need only be relinked with the 
      latest libraries from the distribution to read the new database 
      structure. Calling arguments and user software won't need to change -- 
      unless the user wants to use the new information.

   - Automatic management of asteroid & comet databases ... software merge of
      the two database components is managed in the readers such that both can 
      be simultaneously and transparently accessed.

   - Error checking

 Since the reader software is configurable with respect to what data it returns,
 it should be possible to duplicate the output of legacy readers users may have
 for DASTCOM3, reducing the required changes in user code and more rapidly 
 eliminating the need for DASTCOM3.

 If the provided subroutines are not compatible with user application code, the 
 provided "dxlook" tool might be used to retrieve data in a shell or scripted 
 environment. Relevant details are also provided later for implementation of a 
 reader in other languages. The structure of the database records is summarized
 in the text file './dastcom5/doc/dastcom5.map'.


 Building Software and Library
 -----------------------------

 To build the 'dxlook' command-line application, and the library 'libdxread.a' 
 used to link with user application software, move to the source directory, 
 edit the Makefile to set the desired compiler and link flags, then execute a 
 "make all" command. On a Unix/Linux/Mac platform:

      1. cd dastcom5/src

      2. vi Makefile 

         - The 'LIBDXREAD' and 'EXEDXREAD' variables can be altered to point to 
            a different destination for the library archive and executables, if 
            desired.
 
         - 'BITDXREAD' can be altered as necesary to select which memory model 
            (32 or 64-bit) is required for compatibility with user applications.
            Options are '-m32' or '-m64'
 
         - Example settings for four FORTRAN compilers are shown (FC and FFLAGS
            variables). It may be possible to simply alter comment-markers to 
            select the desired compiler. If not, create new FC and FFLAGS as 
            needed.
 
         - Once Makefile is customized for the local system, it can be renamed
            'makefile' (lower-case). This will prevent it from being 
            over-written during future releases while giving local changes 
            execution precedence over the default distribution 'Makefile' 
            (upper-case). However, keep an eye on "news.txt" in case an 
            additional program or subroutine is added to the distribution, and
            the distribution Makefile is changed to include it.

      3. make all

          This will build 'libdxread.a', created as '../lib/libdxread.a'

          Application program 'dxlook' will be created as '../exe/dxlook',
          along with some other example programs in the same directory

          User applications would then link with 'libdxread.a' after 
          compilation to gain access to the database subroutines.

 The library and executable files can be relocated and renamed, though the rest
 of this discussion assumes the above basic directory structure for the purpose
 of example.


 Running dxlook 
 --------------

 'dxlook' is a stand-alone command-line tool to rapidly examine DASTCOM 
 database records.  It can also be used to examine legacy databases such as 
 DASTCOM3 and DASTCOM4. 

 Type '../exe/dxlook' from the command-line to execute the program, assuming 
 the users' working directory location of './dastcom5/src' is unchanged. 

 Type 'help' at the 'dxlook>' prompt to list commands. 

 To view small-body data, enter a record number or, if the Unix/Linux/Mac
 egrep utility is in your path, an asteroid name or designation.  

 Type 'x' to exit the 'dxlook' program. 

 Program directives are:

  Directive           Meaning 
  ------------------  -----------------------------------------------------
  DB     [path/file]  Open specified database(s); DASTCOM 3-5 are supported
  SUMM                Summarize currently opened databases
  INDEX  [path/file]  Specify optional ASCII index file for look-ups
  FIELDS [list]       Set field codes to display (default= -5, all DASTCOM5)
  PAGER  {path}       Toggle output pager w/executable path (default= none)
  LABEL  [ON/OFF]     Toggle display of field labels (default= ON)
  HELP                Display this list of directives (same as "?")
  {Object}            Display {object};  record #, name, designation, SPK ID
  !{X}                Pass {X} to operating system as a command
  X                   Exit program; same as CLOSE, QUIT, EXIT

 'dxlook' is controlled by an input file, keyboard, or combination of the two. 

 It is recommended that set-up information which doesn't change, such as the 
 paths to the latest database files, be placed in an input file called 
 "./dxlook.inp".

 When 'dxlook' is run, it will first look for a file called 'dxlook.inp' in 
 the run-time directory, load anything there if found, then drop into keyboard 
 input, which would typically be specific user-queries of the database. 

 This is what happened during the test above, using an example directive file 
 in the current working 'src' directory, './dastcom5/src/dxlook.inp'

 Create another customized dxlook directive file './dxlook.inp' in the 
 original top-level unzip directory ('cd ../..') containing the lines between 
 the '---' marks:

 ---
 DB    ./dastcom5/dat/dast5_le.dat
 DB    ./dastcom5/dat/dcom5_le.dat
 INDEX ./dastcom5/dat/dastcom.idx
 PAGER /bin/usr/less
 ---

 "PAGER" can be set as needed to point to a preferred text paging program on 
 the local system.

 Now execute './dastcom5/exe/dxlook'. Type 'summ' at the prompt to summarize 
 the databases now open. Type "1" to look-up the object in record #1, "99942" 
 for the object in record 99942, and so on.

 While the default behavior is to display all defined DASTCOM fields, the 
 FIELDS directive can be used to request 'dxlook' display only select data. 
 For example, the directive ...
  
   FIELDS 14,11,802,807,808,809,804,805,806

 ... requests display of only object name (14), designation (11), calendar
 epoch (802), and orbital elements EC, A, QR, W, OM, IN. 

 See the './dastcom5/src/curdef.inc' file for a list of all available fields
 and their numeric codes. 

 Some macros are also defined: 

   'FIELDS -3' requests legacy DASTCOM3 fields only, 
   'FIELDS -4' requests legacy DASTCOM4 fields only, 
   'FIELDS -5' restores the default, which is "all DASTCOM5 fields"

 Note that numeric data is always displayed first, followed by character data, 
 regardless of the order you specify with FIELDS. However, the order specified 
 within those two categories (numeric and character) is maintained for display.

 If the Unix/Linux/Mac 'egrep' command is in your path (type '!which egrep' in
 'dxlook' to determine if it can be found), and INDEX is defined, it will be 
 possible to look up objects based on their names, SPK IDs, designations, and 
 other aliases, including regular expressions.  

 If there are multiple matches from such a search, select the desired object 
 using the unique DASTCOM record number on the left-most part of each index 
 line in the list of matches.

 If 'egrep' is not available in your path, only record numbers may be used
 to do look-ups.  


 Using the FORTRAN Reader Library
 --------------------------------

 The principle behind usage of the reader library is to provide a list of data 
 to be retrieved during future calls to the reader. Each field in an asteroid or
 comet record has a code number, as shown in file './dastcom5/src/curdef.inc'.

 The code for each desired field should be placed in an array in the order in 
 which they are to be returned. This output-request 'map' is then consulted for
 each database read and used to determine what information is to be returned 
 from the read, and in what order.

 There are two primary subroutine calls for programmatic users, DXINI and 
 DXREAD, with several other special-purpose utility routines. 

 DXINI initializes the reader package; the user specifies the database files
 to use and the list of desired information and order in which to return it.

 DXREAD then retrieves one small-body data record, returning the requested 
 data in the specified order from the set-up call previously made to DXINI.

 See './dastcom5/src/dxread.f' comments for a list of all subroutines, details 
 of their calling arguments and specific discussion of how to use them. 

 Some macro field codes ("-3","-4","-5") are available that emulate the return 
 of the old D3READ subroutine for DASTCOM3 and DASTCOM4. They can be used if 
 necessary for a "quick and dirty" inclusion of DASTCOM5, by those previously 
 using D3READ and wanting to minimize changes to their existing software, but 
 are not recommmended for general use. 

 The obsolete D3READ subroutine they emulate returned different data in the 
 same array slots depending on whether the object was a comet or asteroid.

 If upgrading to DASTCOM5, it would be better to develop a list of every piece
 of information needed for both comets and asteroids and put them in the 
 request list. The information will therefore be returned in the same position 
 whether the record is for an asteroid or a comet.  Fields not defined for a 
 particular type of object will be returned as zero or blank, but the position 
 of the data won't change using this approach.

 The following are only four general examples, most relevant to legacy DASTCOM3
 users transitioning to DASTCOM5. Once the basic framework is in place 
 (example #1), minor variations allow adapting to a variety of other purposes 
 (examples #2-4).

 Example 1: Open DASTCOM5 and return all information in record (nominal case)
 ---------
C
C** Parameterize array dimensions for convenience. Might instead use
C   "INCLUDE dxparms.inc" and package NUMNS and NUMCH parameters for
C   better future-proofing
C
      INTEGER      NUMNS0, NUMCH0
      PARAMETER( NUMNS0= 142) ! Max.# of numeric fields to be retrieved
      PARAMETER( NUMCH0=  14) ! Max.# of character fields to be retrieved
C
C** Declare necessary variables for DXINI
      CHARACTER*1  DBNAM(2)*256 
      INTEGER      IR8ORD(NUMNS0), ICHORD(NUMCH0), NR8, NCH, ISTAT
      LOGICAL      BUF, WARN
C
C** Declare necessary variables for DXREAD
      CHARACTER*1  CHOUT5(NUMCH0)*80, CHOUT3*217, CERRMS*340
      INTEGER      IOBJ, IZONE, LSRC, LERR
      REAL*8       R8OUT(NUMNS0)
C
C** Declare local variables
      INTEGER      I
C
C** Initialization settings
      DATA         DBNAM    / 
     &              './dastcom5/dat/dast5_le.dat',  ! Database #1
     &              './dastcom5/dat/dcom5_le.dat' / ! Database #2
      DATA         IR8ORD   / -5,141*0 / ! Macro specifying all DASTCOM5 fields
      DATA         ICHORD   /     14*0 / ! Requested character fields
      DATA         NR8,NCH  / 142, 14  / ! Max. # of num. & char. fields to get
      DATA         BUF,WARN / .F., .F. / ! Normal values
C
C** Initialize DASTCOM reader package
      CALL DXINI( DBNAM, IR8ORD, NR8, ICHORD, NCH, BUF, WARN, ISTAT )
C
C** Check initialization return status
      IF ( ISTAT .NE. 0 ) THEN
       PRINT *,'Error on DXINI, ISTAT= ', ISTAT
       CALL DXERR( CERRMS, LERR )
       PRINT *,CERRMS(1:LERR)
       STOP
      ELSE
       PRINT *,'Nominal initialization'
      END IF

C** Set IOBJ for the desired objects' logical record & retrieve data
      IOBJ= 4 ! Numbered asteroid Vesta
      CALL DXREAD( IOBJ, IZONE, LSRC, R8OUT, CHOUT5, ISTAT )
C
C** Check return status and display data
      IF ( ISTAT .NE. 0 ) THEN
       PRINT *,'Error on DXREAD(), ISTAT= ', ISTAT
       CALL DXERR( CERRMS, LERR )
       PRINT *,CERRMS(1:LERR)
       STOP
      ELSE
       PRINT *,'IZONE= ',IZONE     ! Database zone of logical record IOBJ
       PRINT *,'LSRC = ',LSRC      ! Length of SRC vector
       DO I= 1, NR8                ! Objects' numeric data, ordered by IR8ORD()
        PRINT *,' R8OUT(',I,')= ',R8OUT(I)
       END DO
       IF ( NCH .GT. 1 ) THEN      ! DASTCOM5 character array
        DO I= 1, NCH               ! Objects' character data, order by ICHORD()
         PRINT *,' CHOUT5(',I,')= ',CHOUT5(I)
        END DO
       ELSE
        PRINT *,' CHOUT3= ',CHOUT3 ! Legacy character block
       END IF
      END IF

 Depending on how elaborate user software for accessing DASTCOM is, users may
 want to formally INCLUDE ./dastcom5/src/dxparms.inc in application programs 
 that use DXREAD(), and use the parameterized NUMNS and NUMCH values to enhance 
 compatibility with any future release. 

 However, if the application will only ever use DASTCOM in one way, and 
 will not be concerned with future data that may or may not be added, this
 generalization is not necessary. Further, if only a limited number of 
 fields are of interest, smaller dimensions can be specified, as per DXREAD 
 and DXINI documentation. 


 Example 2: Open DASTCOM5, return only DASTCOM3 [asteroid] information (subset)
 ---------

 Replace three corresponding lines in nominal example #1 above with:

      DATA         IR8ORD   /            ! Request DASTCOM3 asteroid num. fields
     &              201,801,802,803,804,805,806,807,808,810,
     &              811,809,432,433,401,402,439,431,438, 123*0  /
      DATA         ICHORD   /            ! Request DASTCOM3 asteroid char fields
     &               14,1,13,11,4,6,7,8,3,5*0 /
      DATA         NR8,NCH  /  19, 9   / ! Max. # num. & char. fields to get

 To generalize and support potential cometary data return for this case, one 
 could define alternate IR8ORD and ICHORD arrays for DASTCOM3 comet records, 
 then test whether object is comet or asteroid and re-initialize with DXCLOS 
 and DXINI calls to use the type-specific list.

 Even better (2b): define a single list that includes all DASTCOM3 asteroid and 
 comet fields. This avoids having to test what the object is and reinitialize 
 with a new list. 

      DATA         IR8ORD   /            ! Request DASTCOM3 ast+com num. fields
     &              201,801,802,803,804,805,806,807,808,810,
     &              811,809,432,433,
     &              401,402,439,431,438,       ! asteroid unique data
     &              408,409,403,404,151,118*0/ ! stick comet-unique data at end
      DATA         ICHORD   /            ! Request DASTCOM3 ast+com char fields
     &               14,1,13,11,4,6,7,8,3,
     &               9,10, 2*0 /               ! stick comet unique data at end
      DATA         NR8,NCH  /  24, 11  / ! Max. # num. & char. fields to get



 Example 3: Open DASTCOM5, return only DASTCOM3 data in legacy D3READ order 
 ---------

 Replace three corresponding lines in nominal example #1 above with:

      DATA         IR8ORD   / -3,141*0 / ! Macro specifying only DASTCOM3 fields
       .
      DATA         NR8,NCH  /  18,  1  / ! Max. # num. & char. fields to get
       .
      CALL DXREAD( IOBJ, IZONE, LSRC, R8OUT, CHOUT3, ISTAT )

 

 Example 4: Open legacy DASTCOM3, return data in legacy subroutine D3READ order
 ---------

 Replace corresponding lines in nominal example #1 above with:

      DATA         DBNAM    /
     &              '/home/user/data/DASTCOM3' !Legacy database (change pointer)
     &              '                            ' / 
      DATA         IR8ORD   / -3,141*0 / ! Macro specifying only DASTCOM3 fields
       .
      DATA         NR8,NCH  /  18,  1  / ! Max. # num. & char. fields to get
       .
      CALL DXREAD( IOBJ, IZONE, LSRC, R8OUT, CHOUT3, ISTAT )


PROGRAMMER NOTES
 
 Overview of Database Structure
 ------------------------------

 Unlike legacy DASTCOM3 and DASTCOM4 databases, DASTCOM5 is provided as two
 files; "dast5" contains numbered and unnumbered asteroids, while "dcom5"
 contains comet records.

 Division into two files for asteroids and comets is done because the two types
 of bodies contain different information. Comet records are fewer but contain 
 more potential data fields and so are individually larger than the much more
 numerous asteroid records.

 Since a direct access file has one fixed record size, it is advantageous to use
 two separate files for asteroids and comets, each with a different fixed record
 size mediated in software. This avoids having total physical database size 
 driven by the few comets with larger records and reduces empty storage space as
 the number of cataloged asteroids having records smaller than comets grows.

 As a result, while the total physical size of DASTCOM5 is about twice that of 
 DASTCOM3, it can store about seven times as much data.
 
 Both files contain a header record as physical record #1. These header
 records include structural information that permits finding objects within 
 the file.


 DAST5
 -----

 'dast5' contains two groups of objects: all the numbered asteroids followed 
 by all the unnumbered asteroids.

 'Numbered asteroids' are those sequentially assigned a number by the IAU
 Minor Planet Center, as measurements accumulate and knowledge of their orbit 
 becomes robust. 

 'Unnumbered asteroids' have less well-determined orbits, and are typically 
 more recently discovered. Once enough information is available for an 
 unnumbered asteroid to secure its orbit solution, it is numbered, and 
 relocated in the database to the numbered zone.

 Consequently, only numbered asteroid record numbers are fixed. All other
 record numbers can shift as objects are moved around the database, thus the 
 need for an index to support look-ups.

 Record "byte" maps showing the specific content of header, asteroid, and comet
 records are in the file './dastcom5/doc/dastcom5.map'. 

 The 'dast5' file overall structure map is shown below, where BIAS() is as 
 returned by a call to DXBND: 
 
     physical                                     logical
       record  ----------------------              record
            1  Header (record #1)
               ----------------------  --------------------_
    1-BIAS(1)  Numbered asteroid 1                     1    |
    2-BIAS(1)  Numbered asteroid 2                     2    |
                        .                                   |  Zone 1 
                        .                                   |  (numbered ast.)
    N-BIAS(1)  Numbered asteroid N                     N   _|
                        .                                  _
  N-BIAS(1)+1  Unnumbered asteroid #1  N-BIAS(1)+1+BIAS(2)  |
  N-BIAS(1)+2  Unnumbered asteroid #2  N-BIAS(1)+2+BIAS(2)  |
                        .                                   | Zone 2 
                        .                                   | (unnumbered ast.)
  N-BIAS(1)+X  Unnumbered asteroid #X  N-BIAS(1)+X+BIAS(2) _|


 DCOM5
 -----

 The 'dcom5' overall file structure map is simpler and looks like this:
 
     physical                                     logical
       record  ----------------------              record
            1  Header (record #1)
               ----------------------  --------------------_
            2  Comet record #1                   2+BIAS(3)  |
            3  Comet record #2                   3+BIAS(3)  |
                        .                                   |  Zone 3 (comets)
                        .                                   |
          Y+1  Comet record #Y                 Y+1+BIAS(3) _|


 ADDRESSING
 ----------
 
 Considering the level of database flux, with objects being added and moved
 around daily, the databases and DXREAD library package are organized to use 
 "logical record" addressing to look up objects. This offers some flexibility.

 The relationship between physical database record numbers and the logical 
 record numbers users pass to DXREAD to retrieve those records is ...

             physical_record = logical_record - bias_for_zone

 The bias pointer for each zone is given in the header, while the index file 
 provides each objects' logical record number. The above relationship then 
 permits calculation of the object's physical record, which can then be 
 directly read.

 This makes it possible for users to request IAU numbered object 1 (Ceres), for
 example, by simply passing IOBJ= 1 to DXREAD(), even though the data is not 
 physically located in record 1. The reader software can relate that input to 
 the physical database given the structural information in the header. 

 For unnumbered asteroids and comets, the logical address approach allows data 
 to be relocated as necessary for database maintenance while reducing the 
 apparent movement of records visible to users.


 SCANNING THE DATABASE
 ---------------------

 Header records also report the logical record boundaries for the categories of 
 objects. This information can be used, for example, to perform a sequential
 search of the entire database as follows:

 Call DXINI to initialize, then call DXBND to retrieve the logical record 
 bounds (and bias pointers) for each category of object. Proceed to read each 
 record in sequence, calling DXREAD as IOBJ is incremented positively from 
 BND(4) to BND(1) for numbered asterods, BND(5) to BND(2) for unnumbered 
 asteroids, and from BND(6) to BND(3) for comets.


 INDEX FILE
 ----------

 The index file './dastcom5/dat/dastcom.idx' is part of the "matched set" of 
 databases and corresponds to the 'dast5' and 'dcom5' of the same distribution
 archive. The index is also valid for a legacy DASTCOM3 or DASTCOM4 database 
 corresponding to the same DASTCOM5 files.

 The index is a plain-text ASCII file suitable for pattern matching using one 
 of the 'grep' family of tools under Unix/Linux/Mac.

 An example index entry line (between the '----' markers):

----
751439 2013 LT14               ,3640966,2013 LT14,2010 NY13,3536865,K13L14T,
----

 From left to right:

  1: DASTCOM logical record number -- the value to pass to DXREAD() as IOBJ
      (potentially up to an 8 digit integer, for DASTCOM5)

  2: Space

  3: Name followed by a comma (alternatively, primary designation, if not named)

  4: Primary SPK ID number followed by a comma

  5: Primary designation followed by a comma

  6: Variable length list of historical or alternate designations 
      and SPK IDs (from when it was an unnumbered asteroid, for example) 

  7: Final item on the line is an MPC packed designation
 
 If developing a parsing routine for the index, robustness requires keeping   
 the length of each item flexible and using the space or comma delimiters as 
 bounds. In the future, names may grow longer, record and SPK-ID range may 
 change, etc.


 Alternative Access Methods
 --------------------------

 It is strongly recommended the provided readers be used. However, if this is
 not possible due to a particular language requirement or production 
 environment ...

 1. Note that 'dxlook' (perhaps with directive "LABELS off") could be used in 
     an interpreted script or shell environment to do the database interaction,
     with the output being loaded into script variables.

 2. The databases themselves are not FORTRAN specific, but can be read using 
     any language that can load a byte-pattern into a properly typed variable
     (Perl, etc.). Consult the byte map in './dastcom5/doc/dastcom5.map' to 
     determine the size and arrangement of each variable in header, asteroid, 
     and comet records, and have at it. Comments in the FORTRAN code may help. 

     The database records consist of 1-byte, 2-byte, and 4-byte integers, 
     4-byte and 8-byte floating point reals, and character data-types 
     (1-byte sequences).

     The nominal distribution files 'dast5_le.dat' and 'dcom5_le.dat' are 
     little-endian.  Therefore, for multi-byte patterns, the first byte is 
     smallest. 

     The provided readers automatically translate on big-endian systems, but 
     this would have to be handled by user-developed readers on big-endian 
     systems. However, some of those architectures (generally RISC-based) or
     compilers (Intel FORTRAN, for example) may offer switchable endianness. 


CONTENTS OF DASTCOM5

   The following is excerpted from '../src/curdef.inc', which will always be 
the definitive listing:

C CURDEF.INC (FORTRAN)
C  Add definitions to end of lists here if the database is altered to include 
C  additional parameters, but do not change or delete existing assignments.
C
C Data-fields currently defined and retrievable by calling DXREAD are listed
C below. Not all retrievable fields will be populated with data.
C
C Column "CODE" gives the integer used to request that quantity by its placement
C in the IR8ORD() array passed to DXINI. The requested value will be stored in
C the corresponding position of the R8OUT output array returned by DXREAD.
C
C All numeric data is returned to the user by DXREAD as 8-byte IEEE floating
C point, although the numeric data is physically stored in the database with 
C the different byte-level representations indicated below to reduce physical
C storage space. 
C 
C Symbols in the availability code:
C
C   a = field is defined for asteroid records
C   c = field is defined for comet records
C   3+= field is defined for DASTCOM3 and later database readers
C   4+= field is defined for DASTCOM4 and later database readers
C   5+= field is defined for DASTCOM5 and later database readers
C
C The classical orbital elements are heliocentric ecliptic. The ecliptic is
C specified by EQUNOX, which should now always be '2000', indicating J2000
C ecliptic system (IAU76 constants). 
C
C To transform from J2000 ecliptic to the original equatorial system used for 
C numerical integration (defined by PENAM), use an obliquity angle (epsilon)
C of 84381.448 arcsec. There are slightly different values of epsilon in
C different sets of constants, but this value was used to convert from the 
C precise planetary ephemeris equatorial system of PENAM, and so should be
C used to consistently recover the original numerically integrated equatorial 
C coordinates. If ecliptic coordinates are to be used directly for high 
C precision applications, convert-as-necessary to the intended ecliptic system.
C
C If EQUNOX is instead '1950', an epsilon of 84404.8362512 arc-seconds should
C be used to recover the original integrator equatorial state in the coordinate 
C frame of PENAM.
C
C The au->km conversion value for 'au' distance units should be taken from the 
C planetary ephemeris indicated by PENAM, along with other constants such 
C as mass parameters (GM's). For PENAMs of at least DE430 and later, 1 au is 
C 149597870.700 km, by IAU standard.
C
C Time units are in the coordinate time (CT) scale of general relativity, 
C defined by the independent variable in the barycentric planetary ephemeris 
C equations of motion. This is equivalent to the IAU TDB time-scale.
C
C "1 day" is defined to be 86400 SI seconds.
C
C NOTE #1: non-gravitational parameters A1, A2, A3 are stored in the database
C in units of 10^-8 au/d^2, but returned by this reader package with units of
C au/d^2. Conversion is done in subroutine DXNASGN. This handling is denoted 
C in the list below by the "[s:10^-8 au/day^2]" notation, where "s:" denotes 
C "stored as". For example, if the database store value of A1 is 1.3453, the
C value of A1 is 1.3453*10^-8 au/day^2, and that will be returned by this 
C reader package. A user-developed reader would need to perform the same units
C conversion.
C
C NOTE #2: values are physically stored with minimum appropriate byte-lengths 
C in suitable numerical data-types (i.e., integer, floating point), but are 
C converted to an array of 8-byte floating-point double-precision values when 
C returned to the user by DXREAD(). This can result in additional digits that 
C aren't meaningful; for example, when 4-byte floating-point values (seven
C decimal digits) are stored in the 8-byte return array (fifteen decimal 
C digits), the eight additional digits can be artifacts of floating point 
C representation. Calculations should take this into account and recognize 
C the original REAL*4 database values as having at most seven decimal digits 
C of precision.
C
C NOTE #3: The non-gravitational A1-A3, DT and related model parameters are as
C described in "Cometary Orbit Determination and Nongravitational Forces", 
C D.K. Yeomans, P.W. Chodas, G. Sitarski, S. Szutowicz, M Krolikowska, in 
C Comets II, University of Arizona Press (2004), pp. 137-151.
C
C Physically stored values REAL*8 (8-byte IEEE floating precision numeric)
C
C CODE avail   Label  Definition
C ---- ----- -------  ---------------------------------------------------------
C  801 ac/3+   EPOCH  Time of osc. orbital elements solution, JD (CT,TDB)
C  802 ac/3+  CALEPO  Time of osc. orbital elements solution, YYYYDDMM.ffff
C  803 ac/3+      MA  Mean anomaly at EPOCH, deg (elliptical & hyperbolic cases,
C                       "9.999999E99" if not available)
C  804 ac/3+       W  Argument of periapsis at EPOCH, J2000 ecliptic, deg.
C  805 ac/3+      OM  Longitude of ascending node at EPOCH, J2000 ecliptic,deg.
C  806 ac/3+      IN  Inclination angle at EPOCH wrt J2000 ecliptic, deg.
C  807 ac/3+      EC  Eccentricity at EPOCH
C  808 ac/3+       A  Semi-major axis at EPOCH, au
C  809 ac/3+      QR  Perihelion distance at EPOCH, au
C  810 ac/3+      TP  Perihelion date for QR at EPOCH, JD (CT,TDB)
C  811 ac/3+   TPCAL  Perihelion date for QR at EPOCH, format YYYYMMDD.fff
C  812 ac/5+  TPFRAC  Decimal (fractional) part of TP for extended precision
C  813 ac/4+  SOLDAT  Date orbit solution was performed, JD (CT,TDB)
C
C DERIVED from physically stored REAL*8 (8-byte IEEE floating precision numeric)
C
C CODE avail   Label  Definition
C ---- ----- -------  ---------------------------------------------------------
C  851 ac/3+   ADIST  Aphelion distance at EPOCH, au
C  852 ac/3+     PER  Sidereal orbit period @ EPOCH, Julian yr (365.25 d/Jul.yr)
C  853 ac/3+  ANGMOM  Specific angular momentum at EPOCH, au^2/D
C  854 ac/3+       N  Mean motion, deg/day (elliptical and hyperbolic cases,
C                       "9.999999E99" if not available)
C  855 ac/3+     DAN  Heliocentric distance of ascending node at EPOCH, au
C  856 ac/3+     DDN  Heliocentric distance of descending node at EPOCH, au
C  857 ac/3+       L  Ecliptic longitude of perihelion at EPOCH, deg.
C  858 ac/3+       B  Ecliptic latitude of perihelion at EPOCH, deg.
C
C Physically stored REAL*8 (8-byte IEEE floating precision numeric) 
C Requires multiple storage slots 
C
C CODE avail   Label  Definition
C ---- ----- -------  ---------------------------------------------------------
C  899 ac/4+ SRC(01)  Square root covariance vector. Vector-stored upper-
C                     triangular matrix with order {EC,QR,TP,OM,W,IN,{ESTL}}.
C                     Always reserve enough space (i.e., up to 55 slots) in 
C                     both request and output arrays to hold complete vector.
C
C                     SRC in matrix form (units are days, radians, au):
C
C                            EC    QR    TP    OM    W      IN  {ESTL}
C                        EC  e11   e12   e13   e14  e15    e16
C                        QR        e22   e23   e24  e25    e26
C                        TP              e33   e34  e35    e36
C                        OM                    e44  e45    e46
C                        W                          e55    e56
C                        IN                                e66
C                     {ESTL}
C 
C                     The SRC vector stores the matrix elements from row i and 
C                     column j (i<=j) in vector slot j*(j-1)/2 + i. i.e., upper 
C                     part only, stored by columns: 
C                           
C                              SRC(1,..)= e11,e12,e22,e13,e23 ...
C
C                     To obtain a covariance matrix, bottom-fill an SRC matrix
C                     with zeros to obtain matrix RI. Multiply RI by its
C                     transpose:
C
C                                           T
C                              COV = RI * RI
C
C Physically stored values REAL*4 (4-byte IEEE floating precision numeric)
C
C CODE avail  Label   Definition
C ---- ----- -------  ---------------------------------------------------------
C  401 ac/3+       H  Absolute visual magnitude (IAU H-G system) (99=unknown)
C  402 ac/3+       G  Mag. slope parm. (IAU H-G)(99=unknown & 0.15 not assumed) 
C  403  c/3+      M1  Total absolute magnitude, mag.
C  404  c/3+      M2  Nuclear absolute magnitue, mag.
C  405  c/4+      K1  Total absolute magnitude scaling factor
C  406  c/4+      K2  Nuclear absolute magnitude scaling factor
C  407  c/4+   PHCOF  Phase coefficient for K2= 5
C  408 ac/3+      A1  Non-grav. accel., radial component, [s:10^-8 au/day^2]
C  409 ac/3+      A2  Non-grav. accel., transverse component,[s:10^-8 au/day^2]
C  410 ac/4+      A3  Non-grav. accel., normal component, [s:10^-8 au/day^2]
C  411  c/4+      DT  Non-grav. lag/delay parameter, days
C  412 ac/5+      R0  Non-grav. model constant, normalizing distance, au
C  413 ac/5+     ALN  Non-grav. model constant, normalizing factor
C  414 ac/5+      NM  Non-grav. model constant, exponent m
C  415 ac/5+      NN  Non-grav. model constant, exponent n 
C  416 ac/5+      NK  Non-grav. model constant, exponent k
C  417  c/4+      S0  Center-of-light estimated offset at 1 au, km
C  418  c/5+     TCL  Center-of-light start-time offset, d since "ref.time"
C  419 a /5+     LGK  Surface thermal conductivity log_10(k), (W/m/K)
C  420 ac/5+     RHO  Bulk density, kg/m^3
C  421 ac/5+   AMRAT  Solar pressure model, area/mass ratio, m^2/kg 
C  422  c/5+     AJ1  Jet 1 acceleration, au/d^2 
C  423  c/5+     AJ2  Jet 2 acceleration, au/d^2
C  424  c/5+     ET1  Thrust angle, colatitude of jet 1, deg.
C  425  c/5+     ET2  Thrust angle, colatitude of jet 2, deg.
C  426  c/5+     DTH  Jet model diurnal lag angle, deg. (delta_theta)
C  427 ac/5+     ALF  Spin pole orientation, RA, deg.
C  428 ac/5+     DEL  Spin pole orientation, DEC, deg.
C  429 ac/5+  SPHLM3  Earth gravity sph. harm. model limit, Earth radii
C  430 ac/5+  SPHLM5  Jupiter grav. sph. harm. model limit, Jupiter radii
C  431 ac/3+      RP  Object rotational period, hrs
C  432 ac/3+      GM  Object mass parameter, km^3/s^2
C  433 ac/3+     RAD  Object mean radius, km
C  434 ac/5+  EXTNT1  Triaxial ellipsoid, axis 1/largest equat. extent, km
C  435 ac/5+  EXTNT2  Triaxial ellipsoid, axis 2/smallest equat. extent, km
C  436 ac/5+  EXTNT3  Triaxial ellipsoid, axis 3/polar extent, km
C  437 ac/4+    MOID  Earth MOID at EPOCH time, au; '99' if not computed
C  438 ac/3+  ALBEDO  Geometric visual albedo, 99 if unknown 
C  439 a /3+    BVCI  B-V color index, mag., 99 if unknown
C  440 a /5+    UBCI  U-B color index, mag., 99 if unknown
C  441 a /5+    IRCI  I-R color index, mag., 99 if unknown
C  442 ac/4+    RMSW  RMS of weighted optical residuals, arcsec
C  443 ac/5+    RMSU  RMS of unweighted optical residuals, arcsec
C  444 ac/5+    RMSN  RMS of normalized optical residuals
C  445 ac/5+   RMSNT  RMS of all normalized residuals
C  446 a /5+    RMSH  RMS of abs. visual magnitude (H) residuals, mag.
C  447  c/5+   RMSMT  RMS of MT estimate residuals, mag.
C  448  c/5+   RMSMN  RMS of MN estimate residuals, mag.
C
C Physically stored values INTEGER*4 (4-byte integer numeric)
C
C CODE avail   Label  Definition
C ---- ----- -------  ---------------------------------------------------------
C  201 ac/3+      NO  Logical record-number of this object in DASTCOM
C  202 ac/4+    NOBS  Number of observations of all types used in orbit soln.
C  203 ac/4+ OBSFRST  Start-date of observations used in fit, YYYYMMDD 
C  204 ac/4+ OBSLAST  Stop-date of observations used in fit, YYYYMMDD
C
C Physically stored values INTEGER*1 (1-byte signed integers)
C
C CODE avail   Label  Definition
C ---- ----- -------  ---------------------------------------------------------
C  101 ac/5+  PRELTV  Planet relativity "bit-switch" byte: bits 0-7 are set to
C                      1 if relativity for corresponding planet was computed, 
C                      0 if not. For example, if Earth & Jupiter, FORTRAN(95) 
C                      statement IBITS(PRELTV,J,1) should return 1 when J=2 or 
C                      J=4, but zero for every other J through 7. No provision
C                      for supporting Pluto relativity.
C  102 ac/5+  SPHMX3  Earth grav. model max. degree; 0=point-mass, 2= J2 only,
C                      3= up to J3 zonal, 22= 2x2 field, 33=3x3 field, etc.
C  103 ac/5+  SPHMX5  Jupiter grav. max. deg.; 0=point-mass, 2= J2 only, 
C                      3= up to J3 zonal, 22= 2x2 field, 33=3x3 field, etc.
C  104 ac/5+   JGSEP  Galilean satellites used as sep. perturbers; 0=no 1=yes
C  105 ac/5+  TWOBOD  Two-body orbit model flag; 0=no 1=yes
C  106 ac/5+   NSATS  Number of satellites; 99 if unknown.
C  107 ac/4+   UPARM  Orbit condition code; 99 if not computed
C  108 ac/4+    LSRC  Length of square-root cov. vector SRC (# elements used)
C   
C Physically stored values INTEGER*2 (2-byte integers)
C
C CODE avail   Label  Definition
C ---- ----- -------  ---------------------------------------------------------
C  151  c/3+    IPYR  Perihelion year (i.e., 1976, 2012, 2018, etc.)
C  152 ac/3+    NDEL  Number of radar delay measurements used in orbit soln.
C  153 ac/3+    NDOP  Number of radar Doppler measurements used in orbit soln.
C  154  c/5+  NOBSMT  Number of magnitude measurements used in total mag. soln.
C  155  c/5+  NOBSMN  Number of magnitude measurements used in nuc. mag. soln.
C  156  c/3+  COMNUM  IAU comet number (parsed from DESIG)
C 
C Physically stored character data return in CHOUT argument:
C 'length' is the maximum number of characters in the field
C 
C CODE avail  length  Label  Definition
C ---- ----- ------- ------  ---------------------------------------------------
C  001 ac/3+     4   EQUNOX  Equinox of orbital elements ('1950' or '2000')
C  002 ac/4+     6    PENAM  Planetary ephemeris ID/name
C  003 ac/3+    12    SBNAM  Small-body perturber ephemeris ID/name
C  004 a /3      5   SPTYPT  Tholen spectral type
C  005 a /4+     5   SPTYPS  SMASS-II spectral type
C  006 ac/3+     9     DARC  Data arc span (year-year, OR integer # of days)
C  007 a /3+    41   COMNT1  Asteroid comment line #1
C  008 a /3+    80   COMNT2  Asteroid comment line #2
C  009  c/3+    49   COMNT3  Comet comment line #1
C  010  c/3+    80   COMNT4  Comet comment line #2
C  011 ac/3+    13    DESIG  Object designation
C  012 ac/4+    14     ESTL  Dynamic parameter estimation list. Last symbol set
C                             to '+' if list is too long for field; check 
C                             object record comments field for full list.
C  013 ac/3+    10     IREF  Solution reference/ID/name
C  014 ac/3+    29     NAME  Object name


DASTCOM5 BYTE MAP

 The following is an inclusion of './dastcom5/doc/dastcom5.map', which will
 always be the definitive listing:

                     DASTCOM5 Binary Record Structure

Asteroid Header Record  : 835 bytes | Comet Header Record  : 976 bytes
____________________________________|____________________________________
 Bytes reserved         :  86       |  Bytes reserved      :  82
 Bytes undefined        : 749       |  Bytes undefined     : 894
                                    |
Name      Type   Bytes   Cumulative | Name      Type   Bytes   Cumulative
-------  ------- ------  ---------- | -------  ------- ------  ----------
IBIAS1      I       4        4      | IBIAS2      I       4        4
BEGINP      C  8*3=24       28      | BEGINP      C  8*3=24       28
ENDPT       C  8*3=24       52      | ENDPT       C  8*3=24       52
CALDATE     C      19       71      | CALDATE     C      19       71
JDDATE      D       8       79      | JDDATE      D       8       79
FTYP        C       1       80      | FTYP        C       1       80 
BYTE2A      I       2       82      | BYTE2C      I       2       82     
IBIAS0      I       4       86      |                                
<UNDEFINED> -     749      835      | <UNDEFINED> -     894      976    
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 NOTES on DASTCOM5 headers

  1. ENDPT and BEGINP character strings contained the first and last 
     logical record boundaries for the 3 categories of object: 
     IAU-numbered asteroids, unnumbered asteroids, and comets.

     Because the database is issued in two files (DAST5 and DCOM5),
     comet header BND() is zero-filled for DAST5 and asteroid header 
     BND() is zero-filled for DCOM5. 

       Example:

         DAST5:  BEGINP= '000000010050000100000000'
                          |BND(4)||BND(5)||NOCOMS|

                 ENDPT = '003332730075428200000000'
                          |BND(1)||BND(2)||NOCOMS|

         DCOM5:  BEGINP= '000000000000000000900001'
                          |_NO_ASTEROIDS_||BND(6)|

                 ENDPT = '000000000000000000903900'
                          |_NO_ASTEROIDS_||BND(3)|

  2. CALDATE is character-format creation date: "YYYY-MM-DD_HH:MM:SS"
     JDDATE  is Julian astronomical day creation date (REAL*8 numeric)

  3. IBIAS1 & IBIAS2 and new IBIAS0

      Since the database is issued in two files (DAST5 and DCOM5), there
      is no IBIAS2 field in the DAST5 header and no IBIAS0 or IBIAS1 
      fields in the DCOM5 header.

               physical record = user logical record - IBIAS

  4. For byte #80 of header records (character), 

              'A' denotes DASTCOM3, 
              'C' denotes DASTCOM4, 
              '5' denotes DASTCOM5

     The database type variations will be detected & handled in reader(s).
  
  5. Byte 81-82 of each header is a two-byte integer always set to 26901 
     and necessary for database structure tests.

  6. <UNDEFINED> may be null, character  ' ', or anything

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Asteroid Record : 835 bytes         | Comet Record     : 976 bytes
____________________________________|____________________________________
 Bytes reserved : 835               |   Bytes reserved : 976
 Bytes undefined:   0               |   Bytes undefined:   0
                                    |
Name      Type   Bytes   Cumulative | Name     Type    Bytes   Cumulative 
-------  ------- ------  ---------- | ------  -------  ------  ---------- 
NO          I       4         4     | NO         I       4         4   
NOBS        I       4         8     | NOBS       I       4         8      
OBSFRST     I       4        12     | OBSFRST    I       4        12 
OBSLAST     I       4        16     | OBSLAST    I       4        16   
EPOCH       D       8        24     | EPOCH      D       8        24    
CALEPO      D       8        32     | CALEPO     D       8        32  
MA          D       8        40     | MA         D       8        40  
W           D       8        48     | W          D       8        48  
OM          D       8        56     | OM         D       8        56  
IN          D       8        64     | IN         D       8        64  
EC          D       8        72     | EC         D       8        72  
A           D       8        80     | A          D       8        80  
QR          D       8        88     | QR         D       8        88  
TP          D       8        96     | TP         D       8        96  
TPCAL       D       8       104     | TPCAL      D       8       104  
TPFRAC      D       8       112     | TPFRAC     D       8       112 
SOLDAT      D       8       120     | SOLDAT     D       8       120
SRC(1-45)   D 45*8=360      480     | SRC(1-55)  D 55*8= 440     560
PRELTV      I       1       481     | PRELTV     I       1       561 
SPHMX3      I       1       482     | SPHMX3     I       1       562 
SPHMX5      I       1       483     | SPHMX5     I       1       563 
JGSEP       I       1       484     | JGSEP      I       1       564 
TWOBOD      I       1       485     | TWOBOD     I       1       565 
NSATS       I       1       486     | NSATS      I       1       566 
UPARM       I       1       487     | UPARM      I       1       567
LSRC        I       1       488     | LSRC       I       1       568
                                    | IPYR       I       2       570
NDEL        I       2       490     | NDEL       I       2       572
NDOP        I       2       492     | NDOP       I       2       574
                                    | NOBSMT     I       2       576
                                    | NOBSMN     I       2       578
H           R       4       496     | H          R       4       582
G           R       4       500     | G          R       4       586
                                    | M1 (MT)    R       4       590
                                    | M2 (MN)    R       4       594
                                    | K1 (MTS)   R       4       598
                                    | K2 (MNS)   R       4       602
                                    | PHCOF (MNP)R       4       606
A1          R       4       504     | A1         R       4       610
A2          R       4       508     | A2         R       4       614
A3          R       4       512     | A3         R       4       618
                                    | DT         R       4       622
R0          R       4       516     | R0         R       4       626
ALN         R       4       520     | ALN        R       4       630
NM          R       4       524     | NM         R       4       634
NN          R       4       528     | NN         R       4       638
NK          R       4       532     | NK         R       4       642
                                    | S0         R       4       646
                                    | TCL        R       4       650
LGK         R       4       536     |                               
RHO         R       4       540     | RHO        R       4       654
AMRAT       R       4       544     | AMRAT      R       4       658
                                    | AJ1        R       4       662
                                    | AJ2        R       4       666
                                    | ET1        R       4       670
                                    | ET2        R       4       674
                                    | DTH        R       4       678
ALF         R       4       548     | ALF        R       4       682
DEL         R       4       552     | DEL        R       4       686
SPHLM3      R       4       556     | SPHLM3     R       4       690
SPHLM5      R       4       560     | SPHLM5     R       4       694
RP          R       4       564     | RP         R       4       698
GM          R       4       568     | GM         R       4       702
RAD         R       4       572     | RAD        R       4       706
EXTNT1      R       4       576     | EXTNT1     R       4       710
EXTNT2      R       4       580     | EXTNT2     R       4       714
EXTNT3      R       4       584     | EXTNT3     R       4       718
MOID        R       4       588     | MOID       R       4       722
ALBEDO      R       4       592     | ALBEDO     R       4       726
BVCI        R       4       596     |                               
UBCI        R       4       600     |                               
IRCI        R       4       604     |                               
RMSW        R       4       608     | RMSW       R       4       730
RMSU        R       4       612     | RMSU       R       4       734
RMSN        R       4       616     | RMSN       R       4       738
RMSNT       R       4       620     | RMSNT      R       4       742
RMSH        R       4       624     |                               
                                    | RMSMT      R       4       746
                                    | RMSMN      R       4       750
EQUNOX      C       4       628     | EQUNOX     C       4       754
PENAM       C       6       634     | PENAM      C       6       760
SBNAM       C      12       646     | SBNAM      C      12       772
SPTYPT      C       5       651     |                               
SPTYPS      C       5       656     |                               
DARC        C       9       665     | DARC       C       9       781
COMNT1      C      41       706     |                               
COMNT2      C      80       786     |                               
                                    | COMNT3     C      49       830
                                    | COMNT2     C      80       910
DESIG       C      13       799     | DESIG      C      13       923
ASTEST      C       8       807     | COMEST     C      14       937
IREF        C      10       817     | IREF       C      10       947
ASTNAM      C      18       835     | COMNAM     C      29       976
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  'D' type denotes floating-point DOUBLE PRECISION (8-byte REAL*8)
  'I' type denotes INTEGER of indicated byte-length (1,2,or 4)
  'R' type denotes floating-point REAL*4 (4-byte)
  'C' type denotes CHARACTER data of indicated length
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


CONTACT
 
 If any bugs are found, or problems are encountered compiling or using the 
 FORTRAN code provided, with compiler X on platform Y, let us know: 

 e-mail: Jon.D.Giorgini@jpl.nasa.gov

 If you are able to successfully build and link on a different OS with a
 different compiler, let us know how it went and what you had to do. Credit 
 will be given where credit is due.

 This software is primarily intended for our internal use (and we aren't a
 programming group) so is being made available "as-is". However, some issues 
 might readily be resolved if access to compiler X on platform Y is possible.

 Little or no support is available to assist re-implementation outside JPL in 
 other languages, or resolve problems with significantly altered versions of 
 the distribution. 

COPYRIGHT

Copyright 2013, by the California Institute of Technology. ALL RIGHTS RESERVED. 
United States Government Sponsorship acknowledged. Any commercial use must be 
negotiated with the Office of Technology Transfer at the California Institute 
of Technology.
 
This software may be subject to U.S. export control laws. By accepting this
software, the user agrees to comply with all applicable U.S. export laws and
regulations. User has the responsibility to obtain export licenses, or other
export authority as may be required before exporting such information to
foreign countries or providing access to foreign persons.
