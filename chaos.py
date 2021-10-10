def rossler(x,y,z,h,T):
    xxr = np.matrix([[x],[y],[z]])
    for i in range (1,T):
        x1 =x + h*(-y-z)
        y1 =y + h*(x + 0.1*y)
        z1 =z + h*(0.1 + z*(x-14))
        x=x1
        y=y1
        z=z1
        xn=np.matrix([[x],[y],[z]])
        xxr = np.append(xxr,xn,1)
        i = i+1
    return xxr

def Lorenz(state, t, rho, sigma, beta,h,T):
    xxr = np.matrix([[x],[y],[z]])
    x, y, z = state  # Unpack the state vector
    for i in range (1,T):
        x1, y1, z1= sigma * (y - x), x * (rho - z) - y, x * y - beta * z
        x, y, z = x1, y1, z1
        xn
        xn=np.matrix([[x],[y],[z]])
        xxr = np.append(xxr,xn,1)
        i = i+1
    return xxr

def JacobianLorenz965D(x):
    yo=np.matrix([[1,0.0,0.0,0.0,0.0],[0.0,1,0.0,0.0,0.0],[0.0,0.0,1,0.0,0.0],[0.0,0.0,0.0,1,0.0],
               [0.0,0.0,0.0,0.0,1]])
    y = np.matrix([[-1,x[4,0],0,-x[4,0],x[1,0]-x[3,0]],
                  [x[2,0]-x[4,0],-1,x[0,0],0,-x[0,0]],
                  [-x[1,0],x[3,0]-x[0,0],-1,x[1,0],0],
                  [0,-x[2,0],x[4,0]-x[1,0],-1,x[2,0]],
                 [x[3,0],0,-x[3,0],x[0,0]-x[2,0],-1]])
    ya = yo + 0.001*y
    return (ya)

def Lorenz96(x, dt,N,T):
    # Setting up vector
    xa = x
    d = np.zeros((N,1))
    # Loops over indices (with operations and Python underflow indexing handling edge cases)
    for i in range (1,T):
        for i in range(N):
            d[i] = (x[(i + 1) % N] - x[i - 2]) * x[i - 1] - x[i] + 8
            xn = x + dt*d
            x=xn
        xa = np.append(xa,x,1)   
    return xa

def RungeKutta43d(f,g,p,x0,y0,z0,h,T):
    t = np.linspace(0,T,1000)
    N = len(t)
    # Initialise vectors
    x = np.zeros((N))
    y = np.zeros((N)) 
    z = np.zeros((N))
    # Starting conditions
    x[0] = x0
    y[0] = y0
    z[0] = z0

    # Initialise derivative functions
    #dx = f(x,y,z,dt)                  # dx = x' = dx/dt
    #dy = g(x,y,z,dt)   # dy = y' = dy/dt
    #dz = h(x,y,z,dt)  # dz = z' = dz/dt

    # Initialise K vectors
    kx = np.zeros((4)) # to store K values for x
    ky = np.zeros((4)) # to store K values for y
    kz = np.zeros((4)) # to store K values for z
    b = np.matrix([1, 2, 2, 1])   # RK4 coefficients

    # Iterate, computing each K value in turn, then the i+1 step values
    for i in range(0,N-1):        
        kx[0] = f( x[i], y[i], z[i],t[i])
        ky[0] = g( x[i], y[i], z[i],t[i])
        kz[0] = p( x[i], y[i], z[i],t[i])

        kx[1] = f( x[i] + (h/2)*kx[0], y[i] + (h/2)*ky[0], z[i] + (h/2)*kz[0],t[i] + (h/2))
        ky[1] = g( x[i] + (h/2)*kx[0], y[i] + (h/2)*ky[0], z[i] + (h/2)*kz[0],t[i] + (h/2))
        kz[1] = p( x[i] + (h/2)*kx[0], y[i] + (h/2)*ky[0], z[i] + (h/2)*kz[0],t[i] + (h/2))

        kx[2] = f( x[i] + (h/2)*kx[1], y[i] + (h/2)*ky[1], z[i] + (h/2)*kz[1],t[i] + (h/2))
        ky[2] = g( x[i] + (h/2)*kx[1], y[i] + (h/2)*ky[1], z[i] + (h/2)*kz[1],t[i] + (h/2))
        kz[2] = p( x[i] + (h/2)*kx[1], y[i] + (h/2)*ky[1], z[i] + (h/2)*kz[1],t[i] + (h/2))

        kx[3] = f( x[i] + h*kx[2], y[i] + h*ky[2], z[i] + h*kz[2],t[i] + h)
        ky[3] = g( x[i] + h*kx[2], y[i] + h*ky[2], z[i] + h*kz[2],t[i] + h)
        kz[3] = p( x[i] + h*kx[2], y[i] + h*ky[2], z[i] + h*kz[2],t[i] + h)

        x[i+1] = x[i] + (h/6)*(kx[0]+2*kx[1]+2*kx[2]+kx[3])
        y[i+1] = y[i] + (h/6)*(ky[0]+2*ky[1]+2*ky[2]+ky[3])       
        z[i+1] = z[i] + (h/6)*(kz[0]+2*kz[1]+2*kz[2]+kz[3]) 
        i = i+1
    return(x)