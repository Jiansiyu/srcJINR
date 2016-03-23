// -------------------------------------------------------------------------
// -----                 FairWebScreenshots header file                -----
// -----                Created 11/12/15  by K. Smirnov                -----
// ------------------------------------------------------------------------- 

#ifndef FAIRWEBSCREENSHOTS_H
#define FAIRWEBSCREENSHOTS_H

#include "FairTask.h"
#include "FairEventManager.h"           // for FairEventManager

#include "TString.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <netinet/in.h>

///////////////////////////////////////////////////////////////////////////////
#define BUFFER_SIZE 512
#define MAX_FILE_SIZE 1024*1024
#define MAX_CONNECTIONS 3
///////////////////////////////////////////////////////////////////////////////

struct www_thread_par
{
    int web_port;
    TString outputDir;
    int daemon;
};


class FairWebScreenshots : public FairTask
{    
  public:
    // Standard constructor
    //*@param name        Name of task
    //*@outputDir         Output directory
    //*@param iVerbose    Verbosity level
    FairWebScreenshots(const char* name, char* output_dir, bool isWebServer = false, Int_t iVerbose = 1);

    // Destructor 
    virtual ~FairWebScreenshots();

    // Set verbosity level. For this task and all of the subtasks. 
    void SetVerbose(Int_t iVerbose) { fVerbose = iVerbose; }
    // Executed task 
    virtual void Exec(Option_t* option);
    virtual InitStatus Init();
    virtual void SetParContainers();
    // Action after each event
    virtual void Finish();

    //Set format of saved files
    void SetFormatFiles(int formatFiles); 
    //Set quantity of files
    void SetMultiFiles(int quantityFiles);
    //Set Number port
    void SetPort(int NumberPort);
    
  private:
    // Default constructor
    FairWebScreenshots();
    FairWebScreenshots(const FairWebScreenshots&);
    FairWebScreenshots& operator=(const FairWebScreenshots&);

    static int daemonize();
    static int sendString(const char *message, int socket);
    static void sendHeader(const char *Status_code, char *Content_Type, int TotalSize, int socket);
    static void sendHTML(char *statusCode, char *contentType, char *content, int size, int socket);
    static void sendFile(FILE *fp, int connecting_socket);
    static int scan(char *input, char *output, int start, int max);
    static int checkMime(char *extension, char *mime_type);
    static int getHttpVersion(char *input, char *output);
    static int GetExtension(char *input, char *output, int max);
    static int Content_Lenght(FILE *fp);
    static int handleHttpGET(char *input, TString output_dir, int connecting_socket);
    static int getRequestType(char *input);
    static int receive(int connecting_socket, TString output_dir);
    static int acceptConnection(int currentSocket, TString output_dir);
    static int start(int webPort, TString output_dir);
	static int start_server(void * ptr);

   FairEventManager* fMan;

   // 0 - PNG, 1 -JPG, 2 - both types
   int formatFiles;
   // 0 - one (the same) file event.png; 1 - multiple files event_nnn.png (event1.png and etc.)
   bool quantityFiles;

   TString outputDir;
   int web_port;
   int daemon;

   bool isWebStarted;
   bool isWeb;

 ClassDef(FairWebScreenshots,1);
};

#endif
