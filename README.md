
# Roll Damping Coefficient Estimation (Streamlit App)

This project provides a web-based application to estimate the roll damping coefficient using the IKEA method. It is packaged into a container using Streamlit for an interactive interface.

## Features

- **Streamlit interface**: Easily input data and get results for roll damping coefficient estimation.
- **Dockerized**: The entire application is containerized for ease of deployment.
- **GitHub Integration**: The code is cloned directly from this repository and runs inside the Docker container.

## Getting Started

### Prerequisites

You will need the following installed to run the application locally:

- [Docker](https://docs.docker.com/get-docker/)
  
### Running the Application

1. Clone the repository:

   ```bash
   git clone https://github.com/pciuh/roll-damping-streamlit.git
   cd roll-damping-streamlit
   ```

2. Build the Docker image:

   ```bash
   docker build -t roll-damping-app .
   ```

3. Run the Docker container:

   ```bash
   docker run -p 8501:8501 roll-damping-app
   ```

4. Access the application by opening [http://localhost:8501](http://localhost:8501) in your browser.

### Health Check

The application has a built-in health check, which can be accessed at:

```
http://localhost:8501/_stcore/health
```

This ensures that the Streamlit service is running correctly.

## Structure

- `roll_damp.py`: The main script that runs the roll damping coefficient estimator.
- `requirements.txt`: Lists the Python dependencies required to run the application.
- `Dockerfile`: Defines the container environment.

## Customization

To modify the application:

1. Clone the repository.
2. Make changes to the code as needed.
3. Rebuild the Docker image with the updated code.

```bash
docker build -t roll-damping-app .
```

## License

This project is licensed under the MIT License.

## Contact

For any questions, feel free to reach out via GitHub issues.
