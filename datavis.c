#include <datavis.h>
#include <string.h>
#include <termios.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <pthread.h>
#include "bessel.h"
#include "acs-datagen.h"

volatile sig_atomic_t done = 0;
void sighandler(int sig)
{
    done = 1;
}

/**
 * @brief DataVis data structure.
 * 
 */
data_packet g_datavis_st;
/**
 * @brief Mutex to ensure atomicity of DataVis and ACS variable access.
 * 
 */
pthread_mutex_t datavis_mutex;
/**
 * @brief Condition variable used by ACS to signal to DataVis that data is ready.
 * 
 */
pthread_cond_t datavis_drdy;

typedef struct sockaddr sk_sockaddr;

int main(int argc, char *argv[])
{
    signal(SIGINT, sighandler);

    // init for bessel coefficients
    calculateBessel(bessel_coeff, SH_BUFFER_SIZE, 3, BESSEL_FREQ_CUTOFF);
    // initialize target omega
    z_g_W_target = 1;                       // 1 rad s^-1
    MATVECMUL(g_L_target, MOI, g_W_target); // calculate target angular momentum

    int server_fd, new_socket = -1;
    struct sockaddr_in address;
    int opt = 1;
    int addrlen = sizeof(address);

    // Creating socket file descriptor
    if ((server_fd = socket(AF_INET, SOCK_STREAM, 0)) == 0)
    {
        perror("socket failed");
        pthread_exit(NULL);
    }

    // Forcefully attaching socket to the port 8080
    if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEADDR,
                   &opt, sizeof(opt)))
    {
        perror("setsockopt reuseaddr");
        pthread_exit(NULL);
    }

    if (setsockopt(server_fd, SOL_SOCKET, SO_REUSEPORT,
                   &opt, sizeof(opt)))
    {
        perror("setsockopt reuseport");
        pthread_exit(NULL);
    }
    int flags = fcntl(server_fd, F_GETFL, 0);
    if (flags == -1)
    {
        perror("fcntl");
        pthread_exit(NULL);
    }
    fcntl(server_fd, F_SETFL, flags | O_NONBLOCK);

    address.sin_family = AF_INET;
    address.sin_addr.s_addr = INADDR_ANY;
    address.sin_port = htons(PORT);

    // Forcefully attaching socket to the port 8080
    if (bind(server_fd, (struct sockaddr *)&address,
             sizeof(address)) < 0)
    {
        perror("bind failed");
        pthread_exit(NULL);
    }
    if (listen(server_fd, 3) < 0)
    {
        perror("listen");
        pthread_exit(NULL);
    }
    for (int i = 0; i < 10; i++)
        readSensors();
    memcpy(g_datavis_st.data.start, "FBEGIN", 6);
    memcpy(g_datavis_st.data.end, "FEND", 4);
    while (!done)
    {
        readSensors();
        g_datavis_st.data.tstart = 0;
        g_datavis_st.data.tnow = acs_ct * 100 * 1000;
        // VECTOR_ASSIGN(B, g_datavis_st.data., g_B[mag_index]);
        {
            g_datavis_st.data.x_B = x_g_B[mag_index];
            g_datavis_st.data.y_B = y_g_B[mag_index];
            g_datavis_st.data.z_B = z_g_B[mag_index];
        }
        // VECTOR_ASSIGN(Bt, g_datavis_st.data., g_Bt[mag_index]);
        {
            g_datavis_st.data.x_Bt = x_g_Bt[bdot_index];
            g_datavis_st.data.y_Bt = y_g_Bt[bdot_index];
            g_datavis_st.data.z_Bt = z_g_Bt[bdot_index];
        }
        // VECTOR_ASSIGN(W, g_datavis_st.data., g_W[mag_index]);
        {
            g_datavis_st.data.x_W = x_g_W[omega_index];
            g_datavis_st.data.y_W = y_g_W[omega_index];
            g_datavis_st.data.z_W = z_g_W[omega_index];
        }
        // VECTOR_ASSIGN(S, g_datavis_st.data., g_S[mag_index]);
        {
            g_datavis_st.data.x_S = x_g_S[mag_index];
            g_datavis_st.data.y_S = y_g_S[mag_index];
            g_datavis_st.data.z_S = z_g_S[mag_index];
        }
        char buf[PACK_SIZE + sizeof(char)];
        memcpy(buf + sizeof(char), g_datavis_st.buf, PACK_SIZE);
        *(char *)buf = PACK_SIZE;
        ssize_t sz = send(new_socket, buf, sizeof(buf), MSG_NOSIGNAL);
        if (sz < 0 && !done)
        {
            if ((new_socket = accept(server_fd, (struct sockaddr *)&address,
                                     (socklen_t *)&addrlen)) < 0)
            {
#ifdef SERVER_DEBUG
                perror("accept");
#endif
            }
        }
        usleep(1000000 / 10); // 10 Hz, 100 ms
    }
    close(new_socket);
    close(server_fd);
    return 0;
}